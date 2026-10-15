"""Tk desktop frontend for uploading or pasting ISOVIZ_INPUT documents."""
from __future__ import annotations

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

from ...application.pipeline import (
    PipelineRequest,
    PreparedIsoVIZ,
    dispatch_to_isoviz,
    prepare_isoviz,
)
from ...application.input_selection import document_from_text, load_text_document
from ...application.environment import find_isoviz_launcher
from ...application.input_selection import best_model_root


def format_preparation_report(prepared: PreparedIsoVIZ) -> str:
    """Return a user-readable report without overstating GUI acceptance."""
    report = prepared.patch_report
    validation = prepared.validation
    lines = [
        f"静态检查：{validation['status']}",
        f"匹配模式：{len(report.matched)} / IsoVIZ 共 {validation['mode_count']} 个模式",
        f"振幅来源：{prepared.request.amplitudes.source}",
        f"结构来源：{prepared.request.structure.source}",
    ]
    if prepared.output_path is not None:
        lines.append(f"临时文件：{prepared.output_path}")
    if report.matched:
        lines.extend(["", "已写入的振幅："])
        lines.extend(f"  {amplitude: .5f}  {label}" for label, amplitude in report.matched)
    if report.unmatched_csv:
        lines.extend(["", "CSV 中未使用的模式："])
        lines.extend(f"  {name}" for name in report.unmatched_csv)
    if report.unmatched_isoviz:
        lines.extend(["", "保留原振幅的 IsoVIZ 模式："])
        lines.extend(f"  {label}" for label in report.unmatched_isoviz)
    if validation["errors"]:
        lines.extend(["", "阻止启动的问题："])
        lines.extend(f"  {message}" for message in validation["errors"])
    if validation["warnings"]:
        lines.extend(["", "提醒："])
        lines.extend(f"  {message}" for message in validation["warnings"])
    return "\n".join(lines)


class IsoVIZInputApp:
    """Thin GUI adapter; all scientific processing stays in the pipeline/backend."""

    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("ISOVIZ_INPUT · 振幅导入与启动")
        self.root.geometry("1180x780")
        self.root.minsize(900, 620)
        self.csv_source = tk.StringVar(value="粘贴的 CSV 文本")
        self.structure_source = tk.StringVar(value="粘贴的 .isoviz 文本")
        self.status = tk.StringVar(value="等待输入")
        self._build()

    def _build(self) -> None:
        outer = ttk.Frame(self.root, padding=16)
        outer.grid(row=0, column=0, sticky="nsew")
        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)
        outer.rowconfigure(2, weight=3)
        outer.rowconfigure(5, weight=2)
        outer.columnconfigure(0, weight=1)

        ttk.Label(
            outer,
            text="把振幅输入安全地交给 IsoVIZ",
            font=("Microsoft YaHei UI", 16, "bold"),
        ).grid(row=0, column=0, sticky="w")
        ttk.Label(
            outer,
            text=(
                "CSV 与 .isoviz 都可选择文件或直接粘贴。程序只生成系统临时副本，"
                "不会修改原始输入。"
            ),
        ).grid(row=1, column=0, sticky="ew", pady=(4, 12))

        panes = ttk.Panedwindow(outer, orient=tk.HORIZONTAL)
        panes.grid(row=2, column=0, sticky="nsew")
        csv_frame, self.csv_text = self._input_panel(
            panes,
            title="1. 振幅 CSV",
            source_var=self.csv_source,
            browse_command=self._browse_csv,
            paste_command=lambda: self._paste(self.csv_text, self.csv_source, "粘贴的 CSV 文本"),
            hint="需要 Best Model Parameter（或 amplitude / amp / value）列。",
        )
        structure_frame, self.structure_text = self._input_panel(
            panes,
            title="2. 子群 .isoviz",
            source_var=self.structure_source,
            browse_command=self._browse_structure,
            paste_command=lambda: self._paste(
                self.structure_text,
                self.structure_source,
                "粘贴的 .isoviz 文本",
            ),
            hint="使用与 CSV 相同 IR / OPD / path 的官方子群文件。",
        )
        panes.add(csv_frame, weight=1)
        panes.add(structure_frame, weight=1)

        actions = ttk.Frame(outer)
        actions.grid(row=3, column=0, sticky="ew", pady=(12, 8))
        ttk.Button(
            actions,
            text="检查并生成临时文件",
            command=lambda: self._run(launch=False),
        ).pack(side=tk.LEFT)
        ttk.Button(
            actions,
            text="生成并启动 IsoVIZ",
            command=lambda: self._run(launch=True),
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Label(actions, textvariable=self.status).pack(side=tk.RIGHT)

        launcher = find_isoviz_launcher()
        launcher_text = str(launcher) if launcher else "未发现显式启动器；Windows 将尝试文件关联"
        ttk.Label(outer, text=f"IsoVIZ 启动方式：{launcher_text}").grid(
            row=4,
            column=0,
            sticky="ew",
            pady=(0, 6),
        )
        self.report_text = ScrolledText(outer, height=12, wrap=tk.WORD, state=tk.DISABLED)
        self.report_text.grid(row=5, column=0, sticky="nsew")

    def _input_panel(
        self,
        parent: ttk.Panedwindow,
        *,
        title: str,
        source_var: tk.StringVar,
        browse_command,
        paste_command,
        hint: str,
    ) -> tuple[ttk.LabelFrame, ScrolledText]:
        frame = ttk.LabelFrame(parent, text=title, padding=10)
        frame.rowconfigure(3, weight=1)
        frame.columnconfigure(0, weight=1)
        toolbar = ttk.Frame(frame)
        toolbar.grid(row=0, column=0, sticky="ew")
        ttk.Button(toolbar, text="选择文件…", command=browse_command).pack(side=tk.LEFT)
        ttk.Button(toolbar, text="粘贴剪贴板", command=paste_command).pack(side=tk.LEFT, padx=(6, 0))
        text = ScrolledText(frame, wrap=tk.NONE, undo=True, font=("Consolas", 10))
        ttk.Button(
            toolbar,
            text="清空",
            command=lambda: self._clear(text, source_var, title),
        ).pack(side=tk.LEFT, padx=(6, 0))
        ttk.Label(frame, textvariable=source_var).grid(row=1, column=0, sticky="ew", pady=(8, 2))
        ttk.Label(frame, text=hint).grid(row=2, column=0, sticky="ew", pady=(0, 6))
        text.grid(row=3, column=0, sticky="nsew")
        return frame, text

    def _browse_csv(self) -> None:
        initial = best_model_root()
        kwargs = {"initialdir": str(initial)} if initial.is_dir() else {}
        path = filedialog.askopenfilename(
            title="选择振幅 CSV",
            filetypes=(("CSV", "*.csv"), ("所有文件", "*.*")),
            **kwargs,
        )
        if path:
            self._load_file(self.csv_text, self.csv_source, path, ".csv")

    def _browse_structure(self) -> None:
        path = filedialog.askopenfilename(
            title="选择子群 .isoviz",
            filetypes=(("IsoVIZ structure", "*.isoviz"), ("所有文件", "*.*")),
        )
        if path:
            self._load_file(
                self.structure_text,
                self.structure_source,
                path,
                ".isoviz",
                errors="replace",
            )

    def _load_file(
        self,
        widget: ScrolledText,
        source_var: tk.StringVar,
        path: str,
        suffix: str,
        *,
        errors: str = "strict",
    ) -> None:
        try:
            document = load_text_document(path, expected_suffix=suffix, errors=errors)
        except (OSError, ValueError) as exc:
            messagebox.showerror("无法读取输入", str(exc), parent=self.root)
            return
        widget.delete("1.0", tk.END)
        widget.insert("1.0", document.text)
        source_var.set(document.source)
        self.status.set("文件已载入，可检查或继续编辑")

    def _paste(
        self,
        widget: ScrolledText,
        source_var: tk.StringVar,
        source: str,
    ) -> None:
        try:
            text = self.root.clipboard_get()
        except tk.TclError:
            messagebox.showwarning("剪贴板为空", "没有可粘贴的文本。", parent=self.root)
            return
        widget.delete("1.0", tk.END)
        widget.insert("1.0", text)
        source_var.set(source)
        self.status.set("剪贴板文本已载入")

    @staticmethod
    def _clear(widget: ScrolledText, source_var: tk.StringVar, title: str) -> None:
        widget.delete("1.0", tk.END)
        source_var.set(f"{title}：尚未输入")

    def _run(self, *, launch: bool) -> None:
        self.status.set("正在解析与校验…")
        self.root.update_idletasks()
        try:
            request = PipelineRequest(
                amplitudes=document_from_text(
                    self.csv_text.get("1.0", "end-1c"),
                    source=self.csv_source.get(),
                ),
                structure=document_from_text(
                    self.structure_text.get("1.0", "end-1c"),
                    source=self.structure_source.get(),
                ),
            )
            prepared = prepare_isoviz(request)
            body = format_preparation_report(prepared)
            if launch:
                receipt = dispatch_to_isoviz(prepared)
                body += (
                    "\n\n已向 IsoVIZ 分发："
                    f"{receipt.launcher or 'Windows 文件关联'}\n"
                    "注意：分发成功不等于窗口已打开或数据已被正确识别，请在 IsoVIZ 中核对。"
                )
                self.status.set("已请求启动 IsoVIZ")
            elif prepared.launchable:
                self.status.set("静态检查完成")
            else:
                self.status.set("静态检查失败，未生成可启动文件")
            self._show_report(body)
        except (OSError, RuntimeError, ValueError) as exc:
            self.status.set("处理失败")
            self._show_report(f"处理失败：\n{exc}")
            messagebox.showerror("ISOVIZ_INPUT 处理失败", str(exc), parent=self.root)

    def _show_report(self, text: str) -> None:
        self.report_text.configure(state=tk.NORMAL)
        self.report_text.delete("1.0", tk.END)
        self.report_text.insert("1.0", text)
        self.report_text.configure(state=tk.DISABLED)


def run_gui() -> int:
    """Create the desktop window and enter Tk's event loop."""
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        raise RuntimeError(f"Unable to create the desktop GUI: {exc}") from exc
    IsoVIZInputApp(root)
    root.mainloop()
    return 0
