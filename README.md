# Equation
多元函数, 方程计算库

Depend on ListMatrix
### Equa0.fs
偏导函数和两个选项的 enum。
### Equa1.fs
含有 Class MEqu，
牛顿迭代法求多元非线性方程组。未作病态或者奇异矩阵的处理。
### Equa2.fs
Class MFunc，功能有
+ 梯度下降法求解无约束局部极值 this.Extremum
+ 罚函数法求带不等式约束的极值 this.ExtremumSubjectTo
+ 拉格朗日乘数法求解等式约束的极值 this.ExtremumSubjectTo

Depend on Equa1.fs
