"""features - 中端：按程序功能划分的六个独立功能包。

每个子包只包含实现该部分功能所需的完整代码，通过调用 :mod:`backend` 的通用
底层能力完成自己的任务；以后只修改某个功能时，只需改对应子包。

    input_cif - 解析输入 CIF
    method1   - Method 1（特殊 k 点子群搜索与共用畸变引擎）
    method2   - Method 2（参数 k 的 (3+d)/公度锁定模式计算）
    method3   - Method 3（Inverse Landau / COPL 与 coupled 见证）
    method4   - Method 4（分解、均匀应变与生成回代）
    export    - 四种下载格式 writer 与共享导出合同
"""
