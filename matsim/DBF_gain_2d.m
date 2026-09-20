function G = DBF_gain_2d(fc, theta, phi, array_m, array_n, array_dx, array_dy, theta_0, phi_0, x_in)
%DBF_GAIN_2D 二维矩形阵列数字波束形成(DBF)的归一化阵列增益
%   G = DBF_gain_2d(fc, theta, phi, array_m, array_n, array_dx, array_dy, theta_0, phi_0, x_in)
%
%   输入:
%     fc        - 载波频率 (Hz),用于反算波长 lambda = c/fc
%     theta     - 波束指向俯仰角 (rad),相对 +z 轴,范围 0~pi
%     phi       - 波束指向方位角 (rad),相对 +x 轴,范围 0~2pi
%     array_m   - x 向阵元个数
%     array_n   - y 向阵元个数
%     array_dx  - x 向阵元间距 (m)
%     array_dy  - y 向阵元间距 (m)
%     theta_0   - 信号到达方向俯仰角 (rad),相对 +z 轴
%     phi_0     - 信号到达方向方位角 (rad),相对 +x 轴
%     x_in      - 参考阵元(0,0)收到的复载波信号,1×N
%   输出:
%     G         - 归一化增益 = max(|y_sum|) / (array_m*array_n * max(|x_in|)),范围 [0,1]
%                 (theta,phi) 对准 (theta_0,phi_0) 时 G = 1
%
%   相位约定:theta、theta_0 均从 +z 起算,phi、phi_0 从 +x 起算;
%   信号从 (theta_0,phi_0) 来时,阵元(m,n)相对(0,0)的相位为
%   +2*pi/lambda*(m*dx*sin(theta_0)*cos(phi_0) + n*dy*sin(theta_0)*sin(phi_0));
%   补偿权取到达相位分布的共轭,主瓣指向 (theta,phi)。

c = 3e8;                % 光速 (m/s)
lambda = c/fc;          % 工作波长 (m)

m = 0:array_m-1;        % x 向阵元序号,0 号阵元为相位基准
n = 0:array_n-1;        % y 向阵元序号,0 号阵元为相位基准

% 信号来向在 x/y 方向引入的阵元相位
wx_in = exp(1j*2*pi*array_dx*sin(theta_0)*cos(phi_0)*m/lambda);
wy_in = exp(1j*2*pi*array_dy*sin(theta_0)*sin(phi_0)*n/lambda);

% 波束指向对应的补偿相位
wx = exp(1j*2*pi*array_dx*sin(theta)*cos(phi)*m/lambda);
wy = exp(1j*2*pi*array_dy*sin(theta)*sin(phi)*n/lambda);

% 二维相位矩阵:行对应 y 向阵元,列对应 x 向阵元
w_xy_in = wy_in.' * wx_in;      % array_n × array_m
w_xy    = wy.'    * wx;         % array_n × array_m

% 展开成一维阵元顺序,做时间域扩展
% 第 k 个阵元的接收相位 × 输入信号 -> 第 k 个阵元收到的时域信号
W_in = w_xy_in(:);              % (array_m*array_n) × 1
x_in_array = W_in * x_in;       % (array_m*array_n) × N

% 波束补偿权 = 到达相位分布的共轭,使期望方向同相叠加
W_H = conj(w_xy(:));             % (array_m*array_n) × 1

% 加权合成:残差 e^{j*2*pi/lambda*[...]}
y = W_H .* x_in_array;          % (array_m*array_n) × N
y_sum = sum(y, 1);              % 沿阵元维同相求和,1×N 合成波形

G = max(abs(y_sum)) / (array_m * array_n * max(abs(x_in)));
end
