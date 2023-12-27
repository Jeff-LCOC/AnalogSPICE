# AnalogSPICE

## 项目介绍
AnalogSPICE是一个模拟电路求解器，支持对于包含电阻、直流电压源、直流电流源、MOSFET的模拟电路的DC分析。输入文件为SPICE格式的电路网表，计算结果以命令行输出。

### 使用方法

AnalogSPICE不需要额外的第三方库，因此可以直接使用命令行进行编译

```shell
clang++ -std=c++17 ./source/*.cpp -o analogspice.out
```

运行方法如下：

```shell
./analogspice.out /YourPathToTestSPICE/test.sp
```

### SPICE 语法

注意：节点编号必须为自然数

1. 注释：以`*`开头的行不会被解析
2. 器件模型：暂时只支持MOSFET模型，格式为：`.MODEL <MODEL_ID> VT <VT> MU <μ> COX <COX> LAMBDA <λ> CJ0 <CJ0>`

3. 电阻：`R<tag> <n1> <n2> <value>`
4. 直流电压源： `V<tag> <n1> <n2> DC <value>`
5. 直流电流源：`I<tag> <n1> <n2> DC <value>`
6. MOSFET：`M<tag> <Source> <Drain> <Gate> <type> <Width> <Length> <model_ID>`


