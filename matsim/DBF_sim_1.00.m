clc


% ===== 仿真参数(可更改) =====
fc = 12e6;              % 载频 (Hz),可更改
theta = -1*pi/4;           % 天线主瓣方向(相对法线),可更改
array_num = 8;          % 阵元个数,可更改
theta_0 = pi/2;         % 输入信号到达方向(相对法线),可更改

% ===== 采样参数 =====
fs = 8 * fc;            % 采样率:每个载波周期采 8 个点
K = 128;                % K 个完整载波周期
N = K * (fs/fc);        % 采样点数(整周期数 -> FFT 无泄漏)
t = (0:N-1)/fs;         % 时间轴

% ===== 输入信号:参考阵元(0号)收到的复载波 =====
omega = 2*pi*fc;        % 载波角频率
phi = 0;                % cu shi
x_in = exp(1j*omega*t+phi); % 复指数载波,幅度为 1

% ===== 波束形成:归一化阵列增益 =====
G = DBF_gain(fc, theta, array_num, theta_0, x_in);
fprintf('归一化阵列增益 G = %.6f\n', G);

% ===== 天线方向图:固定主瓣指向 theta,扫描入射角 theta_0 =====
d = 12.8;                         % 阵元间距 (m),与 DBF_gain 内部保持一致
lamda = 3e8/fc;                   % 波长 (m)
theta0_scan_deg = -90:0.1:90;     % 入射角扫描范围(相对法线,deg),0.1° 步进保证零点清晰
theta0_scan = deg2rad(theta0_scan_deg);

% 仿真链路:主瓣指向固定为 theta,入射角逐点扫描,每个角度调一次 DBF_gain
G_pattern = zeros(size(theta0_scan));
for k = 1:numel(theta0_scan)
    G_pattern(k) = DBF_gain(fc, theta, array_num, theta0_scan(k), x_in);
end

% 理论对照:均匀线阵阵因子(Dirichlet 核),u 为补偿后的残余相位步进
u = 2*pi*d*(sin(theta0_scan) - sin(theta))/lamda;
AF_theory = abs(sin(array_num*u/2) ./ (array_num*sin(u/2)));
AF_theory(abs(u) < 1e-9) = 1;    % u→0 时 0/0 取极限 1(主瓣对准点)

% 绘图(对数域:线性域看不出旁瓣;0 值压到 -60 dB 地板)
floor_dB = -60;
G_dB  = 20*log10(max(G_pattern, 10^(floor_dB/20)));
AF_dB = 20*log10(max(AF_theory, 10^(floor_dB/20)));
figure;
plot(theta0_scan_deg, G_dB, 'b-', 'LineWidth', 1.5); hold on;
plot(theta0_scan_deg, AF_dB, 'r--', 'LineWidth', 1.2);
xline(rad2deg(theta), 'k--');    % 主瓣指向参考线
grid on;
ylim([floor_dB, 0.5]);
xlabel('入射角 \theta_0 (^\circ)');
ylabel('归一化增益 (dB)');
title(sprintf('DBF 接收方向图  N=%d, d/\\lambda=%.3f, \\theta=%.0f^\\circ', ...
      array_num, d/lamda, rad2deg(theta)));
legend('仿真 (DBF\_gain)', '理论阵因子', 'Location', 'best');
