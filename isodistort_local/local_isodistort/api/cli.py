"""
命令行 CLI 入口 - isodistort 命令

使用方式：
    isodistort --cif input.cif --subgroup 3 --amplitude 0.1 -o output
"""
import argparse
import sys
from pathlib import Path

from .core_api import IsoDistort


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="isodistort",
        description="本地化 ISODISTORT - 离线晶体畸变分析工具",
    )

    # 输入
    parser.add_argument("--cif", "-i", type=str, required=True,
                        help="输入 CIF 文件路径")
    parser.add_argument("--distortion-type", "-t", type=str,
                        default="displacement",
                        choices=["displacement", "order", "strain", "magnetic"],
                        help="畸变类型 (默认: displacement)")

    # 子群选择
    parser.add_argument("--list-subgroups", action="store_true",
                        help="仅列出所有子群，不生成畸变")
    parser.add_argument("--subgroup", "-s", type=int, default=None,
                        help="子群序号（来自 --list-subgroups）")

    # 畸变参数
    parser.add_argument("--irrep", type=str, default=None,
                        help="指定不可约表示（如 GM4-）")
    parser.add_argument("--amplitude", "-a", type=float, default=1.0,
                        help="畸变幅度 (默认: 1.0)")
    parser.add_argument("--supercell", nargs=3, type=int, default=[1, 1, 1],
                        metavar=("A", "B", "C"),
                        help="超胞大小，三个整数 (默认: 1 1 1)")

    # 输出
    parser.add_argument("--output", "-o", type=str, default="distorted",
                        help="输出文件名（不含后缀，默认: distorted）")
    parser.add_argument("--format", "-f", type=str, action="append",
                        default=["cif"],
                        choices=["cif", "poscar", "xyz"],
                        help="导出格式，可多次指定 (默认: cif)")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    iso = IsoDistort()

    # 1. 加载结构
    print("=" * 50)
    print("步骤 1: 加载结构并识别对称性")
    print("=" * 50)
    iso.load_structure(args.cif)

    # 2. 列出子群
    print("\n" + "=" * 50)
    print("步骤 2: 枚举各向同性子群")
    print("=" * 50)
    iso.list_subgroups(args.distortion_type)

    if args.list_subgroups:
        print("\n仅列出子群，退出。使用 --subgroup N 选择子群并生成畸变。")
        sys.exit(0)

    if args.subgroup is None:
        print("\n错误: 请使用 --subgroup N 指定子群序号")
        print("       或使用 --list-subgroups 查看所有子群")
        sys.exit(1)

    # 3. 选择路径并计算模式
    print("\n" + "=" * 50)
    print("步骤 3: 选择相变路径，计算畸变模式")
    print("=" * 50)
    iso.select_path(args.subgroup, args.distortion_type)

    # 4. 生成畸变
    print("\n" + "=" * 50)
    print("步骤 4: 生成畸变结构")
    print("=" * 50)
    iso.generate_distortion(
        irrep_label=args.irrep,
        amplitude=args.amplitude,
        supercell=args.supercell,
    )

    # 5. 导出
    print("\n" + "=" * 50)
    print("步骤 5: 导出结构文件")
    print("=" * 50)
    iso.export(args.output, formats=args.format)

    print("\n✅ 完成！")


if __name__ == "__main__":
    main()
