# MATLAB Engine for Python - 快速上手示例
import matlab.engine

print("启动 MATLAB...")
eng = matlab.engine.start_matlab()
print(f"MATLAB 已连接，版本: {eng.version()}")

eng.cd(r'd:\projects\radar_cp\matsim', nargout=0)

# === 基础用法（注意：传 double 类型）===

# 1. 调用 MATLAB 函数（参数用 float，不是 int）
x = eng.linspace(0.0, 10.0, 5.0)
print("linspace 结果:", x)

# 2. 直接执行 MATLAB 命令
eng.eval("y = [1,2,3].^2", nargout=0)
y = eng.workspace['y']
print("MATLAB 变量 y =", y)

# 3. 画图
eng.eval("figure; plot([1,2,3,4], [1,4,9,16]); title('Test Plot')", nargout=0)
print("图已在远程窗口打开")

# 4. 查看工作区
eng.eval("whos", nargout=0)

eng.quit()
print("MATLAB 会话已关闭")
