"""
测试：配置加载与环境设置
"""
import os

from isocore.utils import get_config


def test_config_load():
    cfg = get_config()
    # isobyu 二进制路径可解析（存在性由部署保证；WSL 下为 Windows 侧路径）
    assert cfg.iso_bin.name == "iso"
    assert cfg.findsym_bin.name == "findsym"
    assert cfg.temp_dir.exists(), "临时目录已创建"
    assert cfg.output_dir.exists(), "输出目录已创建"
    assert "ISODATA" in os.environ, "ISODATA 环境变量已设置"
    print("配置加载测试通过")
    print(f"   ISODATA = {os.environ['ISODATA']}")


if __name__ == "__main__":
    test_config_load()
