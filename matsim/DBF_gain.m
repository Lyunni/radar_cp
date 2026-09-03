function G = DBF_gain(fc, theta, array_num, theta_0, x_in)
%DBF_GAIN 数字波束形成(DBF)的归一化阵列增益
%   G = DBF_gain(fc, theta, array_num, theta_0, x_in)
%
%   输入:
%     fc        - 载波频率(Hz),用于反算波长 lamda = c/fc
%     theta     - 波束对准方向(rad),相对阵列法线
%     array_num - 阵元个数(>=2)
%     theta_0   - 信号到达方向(rad),相对阵列法线
%     x_in      - 参考阵元(0号)收到的复载波信号,1×N
%   输出:
%     G         - 归一化增益 = max(|y_sum|) / (array_num * max(|x_in|)),范围 [0,1]
%                 theta 对准 theta_0 时 G = 1;失配时按均匀线阵阵因子衰减
%
%   相位约定:theta、theta_0 均从法线起算;信号从 theta_0 来时,阵元 n 相对
%   0 号阵元相位超前 +2*pi*n*d*sin(theta_0)/lamda;补偿权取到达相位分布的
%   共轭(负指数),把超前"减掉",主瓣指向 theta。

c = 3e8;            % 光速 (m/s)
lamda = c/fc;       % 波长 (m)
d = 12.8;           % 阵元间距 (m):固定几何参数,不随 fc 变
                    % fc=12MHz 时 d/lamda≈0.512 略超半波长,扫方向图时注意栅瓣

n = 0:array_num-1;                                          % 阵元序号,0号阵元为相位基准
Wn    = exp(1j*2*pi*n*d*sin(theta)/lamda);                  % 从 theta 来的到达相位分布 1×array_num
Wn_in = exp(1j*2*pi*n*d*sin(theta_0)/lamda);                % 从 theta_0 来的到达相位分布 1×array_num

x_in_array = Wn_in.' * x_in;        % array_num×N:第 n 行 = 第 n 阵元收到的信号
Wn_H = Wn';                         % 相位补偿权 = 到达相位分布的共轭
y = Wn_H .* x_in_array;             % 补偿后残差 e^{j*2*pi*n*d*(sin(theta_0)-sin(theta))/lamda}
y_sum = sum(y, 1);                  % 沿阵元维同相求和,1×N 合成波形

G = max(abs(y_sum)) / (array_num * max(abs(x_in)));   % 归一化增益:对准时为 1
end
