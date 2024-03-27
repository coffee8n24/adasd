clear;clc;
close all
%% 一些参数
Na_ue = 1;         % UE天线数
Rxarray = [16,16];
Na_bs = 256;       % BS天线数
B_bit = 3;      % 移相器比特
M_beam = 4;     % 波束数目
K_stream = 4;
N_rf = 64;
% 射频链路与天线的映射
PS_mat = zeros(Na_bs,N_rf);
for nrf = 1:N_rf
    nrf_2d = [floor((nrf-1)/(Rxarray(2)))+1,(nrf-1)-floor((nrf-1)/(Rxarray(2)))*Rxarray(2)+1];
    PS_mat(nrf_2d(2)+(0:Rxarray(1)-1)*Rxarray(2),nrf) = ones(16,1);
end
Nu = 16;        % 用户个数
% K_stream = 4;  % 用户流数
M = 4;          % 调制进制QPSK
N_iter = 1e1;       % 迭代次数
SCspacing = 30e3;   % 子载波间隔
N_rb = 273;
Nsc = 12*N_rb;       % 子载波数目
fc = 7e9;           % 载波频率
% 能量计算&网格
SNR = (-20:5:20);        % 信噪比(dB)
Es = 1;                 % 符号平均功率归一
% S/N=SNR,I/N=INR
En = Es./10.^(SNR/10);   % 噪声功率
Len_sim = length(En);  % 遍历次数
% 预分配信道内存
% 频域信道Hf：BS天线数*UT天线数*子载波数*用户数
Hf_ul = zeros(Na_bs,Na_ue,Nsc,Nu);
Hf_dl = zeros(Nu,Na_bs,Nsc);
Gf_dl = zeros(Nu,Na_bs,Nsc);
H_rb = zeros(Nu,Na_bs,N_rb);
G_rb = zeros(Nu,Na_bs,N_rb);
% 预分配各用户传输速率
Rate_ue = zeros(Nu,Len_sim,N_iter);
SE_ue = zeros(Nu,Len_sim,N_iter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 主循环
for niter = 1:N_iter
    %% 生成信道
    disp(['第',num2str(niter),'次信道实现'])
    [~,UToffset] = UE_gene(Nu);
    PL_dB = random('Uniform',100,120,Nu,1);
    PL_UT = 10.^(PL_dB/10);
    UToffset(:,1) = random('Uniform',-180-60,-180+60,Nu,1); % 用户角度偏置
%     SNR_UT = Pt_dB-PL_dB-10*log10(En);
    for nu = 1:Nu
        [Hf_ul(:,:,:,nu),~,Ph] = cdlchannel_generate_URA(Rxarray,Na_ue,1,fc,SCspacing,Nsc,UToffset(nu,:));
        % 分布归一化，使每个天线对的能量为1
        for u = 1:Na_bs
        for s = 1:Na_ue
            Pus = norm(squeeze(Hf_ul(u,s,:,nu)),'fro')^2;
            Hf_ul(u,s,:,nu) = Hf_ul(u,s,:,nu)/sqrt(Pus);
        end
        end
        Hf_ul(:,:,:,nu) = Hf_ul(:,:,:,nu)*sqrt(Nsc);
%         Hf_ul(:,:,:,nu) = Hf_ul(:,:,:,nu)/sqrt(Na_bs);
    %     % 均值归一化，使天线对的平均能量为1
    %     Hf(:,:,:,:,nu) = Hf(:,:,:,:,nu)/sqrt(Ph);
    %     squeeze(sum(sum(sum(sum(abs(Hf).^2,1),2),3),4))/N_iter
    end
    P_check = squeeze(sum(sum(sum(abs(Hf_ul).^2,1),2),3))
    % 转置，信道互惠获得下行信道
    for nu = 1:Nu
        for nsc = 1:Nsc
            Hf_dl(nu,:,nsc) = transpose(squeeze(Hf_ul(:,:,nsc,nu)));
        end
        Gf_dl(nu,:,:) = Hf_dl(nu,:,:)/sqrt(PL_UT(nu));
    end
    % 均值获得资源块的信道信息
    for nu = 1:Nu
        for nrb = 1:N_rb
            H_rb(nu,:,nrb) = sum(squeeze(Hf_dl(nu,:,(1:12)+(nrb-1)*12)),2)/12;
        end
        G_rb(nu,:,:) = H_rb(nu,:,:)/sqrt(PL_UT(nu));
    end
%     squeeze(sum(sum(abs(H_rb).^2,3),2))/N_rb*Nsc
    %% 基于信道相关性的用户分组
    % 初始化主要用户
    % 计算信道能量
    P_u = zeros(Nu,1);
    for nu = 1:Nu
        P_u(nu) = norm(squeeze(Hf_dl(nu,:,:)),'fro')^2;
    end
    % 信道能量最大的用户,设置为第一个组的主要用户
    U_center = zeros(1,M_beam);
    [~,U_center(1)] = max(P_u);
    % 求各用户间的相关值
    corr_user = zeros(Nu,Nu);
    for nu1 = 1:Nu
    for nu2 = nu1:Nu
        corr_user(nu1,nu2) = abs(trace(squeeze(Hf_dl(nu1,:,:)).'*conj(squeeze(Hf_dl(nu2,:,:)))))/(norm(squeeze(Hf_dl(nu1,:,:)),'fro')*norm(squeeze(Hf_dl(nu2,:,:)),'fro'));
    end
    end
    corr_user = corr_user'+corr_user-diag(diag(corr_user));
    % 选取与第一个主要用户相关值最小的M_beam-1个用户，作为剩下几个组的主要用户
    [~,corr_order_init] = sort(corr_user(U_center(1),:),'ascend');
    U_center(2:M_beam) = corr_order_init(1:M_beam-1);
    % 进行用户分组与主要用户更新
    flag_end = 1;
    count_group = 0;
    group_index = zeros(Nu,1);
    while(flag_end)
        U_center_cur = U_center;
        group_index(U_center_cur) = (1:M_beam);
        % 除去目标用户外的剩余用户
        user_left = (1:Nu);user_left(U_center_cur) = [];
        % 根据信道相关性最大化，对剩余用户进行分组
        for n_i = 1:length(user_left)
            [~,group_index(user_left(n_i))] = max(corr_user(user_left(n_i),U_center));
        end
        % 更新主要用户
        for n_beam = 1:M_beam
            user_group_temp = find(group_index==n_beam);
            corr_other = sum(corr_user(group_index==n_beam,group_index~=n_beam),2);
            [~,index_temp] = min(corr_other);
            U_center(n_beam) = user_group_temp(index_temp);
        end
        if(U_center == U_center_cur)
            flag_end = 0;
        end
        count_group = count_group + 1;
    end
    %% 子阵列的模拟权
    % 用户组拼接的等效信道
    R_group = zeros(Na_bs,Na_bs,M_beam);    % 子阵列上的组内用户所有载波的协方差
    vec_Wrf = zeros(Na_bs,M_beam);        % 模拟权的向量化
    vec_ang_Wrf = zeros(Na_bs,M_beam);    % 向量化模拟权的相位角
    vec_Wrf_bit_beam = zeros(Na_bs,M_beam);
    for n_beam = 1:M_beam
        % 计算组内用户所有载波的协方差
        user_index = find(group_index==n_beam);
        for nu = 1:length(user_index)
            R_group(:,:,n_beam) = R_group(:,:,n_beam) + squeeze(Hf_dl(user_index(nu),:,:))*squeeze(Hf_dl(user_index(nu),:,:))'/length(user_index);
        end
        % 选取协方差的最大特征值对应的特征向量，归一化后作为模拟权
        [V_group,D_group] = eig(R_group(:,:,n_beam));
        [~,m_inx] = max(diag(D_group));
        vec_Wrf(:,n_beam) = V_group(:,m_inx)./abs(V_group(:,m_inx));
        % 通过最小距离量化实现b比特量化
        vec_ang_Wrf(:,n_beam) = angle(V_group(:,m_inx));
        for nt = 1:Na_bs
            [~,temp_bit] = min(abs(vec_ang_Wrf(nt,n_beam)+pi-2*pi/2^B_bit*(0:2^B_bit-1)));
            vec_Wrf_bit_beam(nt,n_beam) = exp(1i*2*pi*temp_bit/2^B_bit);
        end
    end
    W_analog = zeros(Na_bs,N_rf);
    % 从向量化的模拟权还原为分组结构的模拟权
    for nrf = 1:N_rf
        nrf_2d = [floor((nrf-1)/(Rxarray(2)))+1,(nrf-1)-floor((nrf-1)/(Rxarray(2)))*Rxarray(2)+1];
        W_analog(PS_mat(:,nrf)==1,nrf) = vec_Wrf_bit_beam(PS_mat(:,nrf)==1,nrf_2d(1));
    end
    %% RB资源分配
    RB_U_mat = zeros(Nu,N_rb);
    Hf_eff = zeros(Nu,N_rf,N_rb);
    P_eff = zeros(Nu,N_rb);
    for nrb = 1:N_rb
        Hf_eff(:,:,nrb) = H_rb(:,:,nrb)*W_analog;
        P_eff(:,nrb) = sum(abs(Hf_eff(:,:,nrb)).^2,2);      % 求第nrb个资源块上每个用户的等效信道能量
    end
    [~,inx_sort] = sort(P_eff,1,"descend");    % [B,IX]=sort(A,1,'descend')让A矩阵的每一列进行降序排序，B为降序后的输出矩阵，IX为B在A的索引，即B=A(IX)
    RB_allocation = inx_sort(1:K_stream,:);    % 每个资源块调度K_stream个用户，RB_allocation是K_stream×N_rb的，其记录了每个资源块上被调度的用户索引
    beta_user_rb = zeros(Nu,N_rb);             % 资源块调度变量矩阵
    for nrb = 1:N_rb
        beta_user_rb(RB_allocation(:,nrb),nrb) = ones(K_stream,1);      % 将被调度的用户赋值为1
    end
    Num_RB_ue = sum(beta_user_rb,2);    % 记录每个用户占据的资源块个数
    beta_user_sc = kron(beta_user_rb,ones(1,12));       % 将资源块调度变量矩阵映射成子载波调度变量矩阵
    W_digtal = zeros(N_rf,Nu,Nsc);
    %% ZF预编码
    for nsc = 1:Nsc
        nrb = floor((nsc-1)/12)+1;
        H_eff_rb = H_rb(RB_allocation(:,nrb),:,nrb)*W_analog;       % 等效信道        
        W_digtal(:,RB_allocation(:,nrb),nsc) = H_eff_rb'/(H_eff_rb*H_eff_rb');   % ZF预编码
        W_digtal(:,RB_allocation(:,nrb),nsc) = W_digtal(:,RB_allocation(:,nrb),nsc)/norm(W_analog*W_digtal(:,RB_allocation(:,nrb),nsc),'fro')*sqrt(K_stream);   % ZF预编码
    end
    %% 遍历信噪比
    parfor nsnr = 1:Len_sim       % ？
        %% 功率分配
        % 混合权功率计算
        G_hybrid = zeros(1,Nsc*K_stream);
        for nsc = 1:Nsc
            n_rb = floor((nsc-1)/12)+1;
            G_hybrid((1:K_stream)+K_stream*(nsc-1)) = sum(abs(Hf_dl(RB_allocation(:,n_rb),:,nsc)*W_analog*W_digtal(:,RB_allocation(:,n_rb),nsc)).^2,1);
        end
        % 注水法
        P_allocation = waterfill(Es,En(nsnr)/Nsc./G_hybrid);
%         % 均匀功率分配
%         P_allocation = ones(1,Nsc*K_stream)*Es/(Nsc*K_stream);
        P_all_mat = zeros(Nu,Nsc);
        for nsc = 1:Nsc
            n_rb = floor((nsc-1)/12)+1;
            P_all_mat(RB_allocation(:,n_rb),nsc) = P_allocation((1:K_stream)+K_stream*(nsc-1))';
        end
        P_nk = zeros(Nu,Nu,Nsc);
        %% 接收信号能量
        for nsc = 1:Nsc
        for nu = 1:Nu
            P_nk(nu,:,nsc) = abs(Hf_dl(nu,:,nsc)*W_analog*W_digtal(:,:,nsc)).^2*P_all_mat(nu,nsc);
        end
        end
        sinr_u = zeros(Nu,Nsc);
        for nu = 1:Nu
            user_inter = (1:Nu);user_inter(nu) = [];
            for nsc = 1:Nsc
                sinr_u(nu,nsc) = squeeze(P_nk(nu,nu,nsc))/(sum(P_nk(nu,user_inter,nsc))+En(nsnr)/Nsc);
            end
%             if(Num_RB_ue(nu)~=0)
%                 SE_ue(nu,nsnr,niter) = sum(log2(1+sinr_u(nu,:)))/(Num_RB_ue(nu)*12);
%             end
            SE_ue(nu,nsnr,niter) = sum(log2(1+sinr_u(nu,:)))/Nsc;
            Rate_ue(nu,nsnr,niter) = SCspacing*sum(log2(1+sinr_u(nu,:)));
        end
    end % end of 信噪比遍历
end % end of 主循环
%% 结果处理&绘图
% 各用户速率
Rate_ue_avg = sum(Rate_ue,3)/N_iter;
% 合速率
Rate_sum = sum(Rate_ue_avg,1);
% 各用户频谱效率
SE_ue_avg = sum(SE_ue,3)/N_iter;
% 合谱效
SE_sum = sum(SE_ue_avg,1);
% 绘图
unit_rate = 1e6;
figure(1);
plot(SNR,Rate_sum/unit_rate,'Color','#000000','LineStyle','-','Marker','.','MarkerSize',12,'LineWidth',1.5);
cell_legend = {'合速率'};
hold on
for nu = 1:Nu
    plot(SNR,Rate_ue_avg(nu,:)/unit_rate);
    cell_legend = {cell_legend{:},['用户',num2str(nu)]};
end
xlabel('信噪比 (dB)')
ylabel('速率 (\times 10^6 bps)')
legend(cell_legend)
figure(2);
plot(SNR,SE_sum,'Color','#000000','LineStyle','-','Marker','.','MarkerSize',12,'LineWidth',1.5);
cell_legend = {'合谱效'};
hold on
for nu = 1:Nu
    plot(SNR,SE_ue_avg(nu,:));
    cell_legend = {cell_legend{:},['用户',num2str(nu)]};
end
xlabel('信噪比 (dB)')
ylabel('频谱效率 (bps/Hz)')
legend(cell_legend)