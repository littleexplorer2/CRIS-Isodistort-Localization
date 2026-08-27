"""
isodistort - 晶体畸变分析工具

基于 ISOTROPY 套件（iso + findsym）的离线封装，
实现 ISODISTORT 核心的结构畸变分析能力。

模块分层：
    backend     - iso / findsym 二进制封装（纯群论计算）
    structure   - 晶体结构层：CIF 读写、对称处理、坐标变换
    distortion  - 畸变业务层：相变路径、模式映射、畸变生成
    io          - 结果输出层：结构文件导出
    api         - 对外接口：Python API
    i18n        - English UI strings
    utils       - 工具层：配置、异常、解析工具
    superspace  - (3+d) 超空间内核（nmod = d；IT-C）
"""

__version__ = "0.3.0"
