PAGE = r"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CRIS REA · Local crystallography</title>
<style nonce="__TOKEN__">
:root{font-family:Inter,Segoe UI,Arial,sans-serif;color:#1f3044;background:#f2f5f9}
*{box-sizing:border-box}body{margin:0}header{background:#10283f;color:white;padding:25px max(24px,calc((100vw - 1180px)/2));display:flex;justify-content:space-between;gap:20px}
h1{font-size:24px;letter-spacing:.07em;margin:0 0 6px}header p{margin:0;color:#bed0df;font-size:13px}.local{font-size:12px;background:#244159;border-radius:30px;padding:8px 13px;align-self:center;white-space:nowrap}
main{max-width:1180px;margin:30px auto;padding:0 24px}.target{background:#e2ecf3;border-left:4px solid #338b9b;padding:18px 22px;margin-bottom:25px;border-radius:3px}
.eyebrow{font-size:11px;text-transform:uppercase;letter-spacing:.12em;color:#517183;font-weight:700}
.target h2{font-size:18px;margin:7px 0}.target p{font-size:13px;line-height:1.6;margin:0}
.grid{display:grid;grid-template-columns:minmax(340px,1fr) minmax(340px,1fr);gap:23px}.panel{background:white;border:1px solid #dce4ed;border-radius:9px;padding:23px;box-shadow:0 3px 9px #17324c05}
h2{font-size:17px;margin:0 0 20px}label{display:block;font-size:12px;font-weight:650;margin:16px 0 7px}select,input[type=number],input[type=text],textarea{width:100%;border:1px solid #c8d3de;border-radius:5px;background:#fcfdff;padding:11px;color:#203246}
textarea{height:285px;resize:vertical;font:12px/1.65 Consolas,monospace}small,.hint{font-size:12px;color:#657588;line-height:1.6}.row{display:flex;gap:15px;align-items:end}.row>div{flex:1}input[type=file]{font-size:12px;width:100%}
button{background:#147b8d;border:none;color:white;border-radius:5px;padding:12px 19px;cursor:pointer;font-weight:650;font-size:13px}button:disabled{opacity:.55;cursor:wait}.actions{display:flex;align-items:center;gap:15px;margin-top:20px}
.status{padding:12px 14px;background:#f1f5f9;border-radius:5px;font-size:12px;overflow-wrap:anywhere;line-height:1.6;white-space:pre-line}.status.busy{background:#fff4d7}.status.error{background:#fce9e6;color:#8b352b}
pre{white-space:pre-wrap;overflow-wrap:anywhere;max-height:360px;overflow:auto;background:#10283f;color:#dfedf7;padding:15px;border-radius:5px;font:12px/1.6 Consolas,monospace;min-height:80px}
.output-label{font-size:12px;margin:19px 0 7px;font-weight:650}details{margin-top:20px}summary{cursor:pointer;font-size:12px;font-weight:650}#evidence{background:#f3f6fa;color:#34495c}
footer{margin:24px 0 35px;font-size:12px;color:#607489;line-height:1.7}code{font-family:Consolas,monospace}#input-source{display:block;margin-top:7px}.timeout-field{max-width:115px}.query-builder{border-bottom:1px solid #dce4ed;padding-bottom:18px;margin:18px 0}.query-builder button{margin-top:14px;background:#45677f}
.mode-grid{margin-bottom:28px}.downloads{display:flex;flex-wrap:wrap;gap:8px;margin-top:15px}.downloads button{background:#45677f;padding:9px 12px}.mode-input{height:185px}.section-heading{margin:30px 0 14px}.basis-note{border-left:3px solid #c89c38;padding-left:12px}
@media(max-width:800px){.grid{grid-template-columns:1fr}header{padding:22px 24px}.local{display:none}}
</style>
</head><body>
<header><div><h1>CRIS REA</h1><p>Official ISOTROPY programs · local execution · traceable runs</p></div><span class="local">● &nbsp; Loopback only</span></header>
<main>
<section class="target"><div class="eyebrow">Research priority</div>
<h2>Parent structure + k / IR / OPD → subgroups and distortion modes</h2>
<p>上传独立母相 CIF，指定 IR / OPD 与必要的 k 参数，计算子群允许的位移和应变子空间。下方保留官方程序原始输入控制台。</p>
<p class="basis-note"><strong>适用范围：</strong>普通三维（ordinary 3D）、公度（commensurate）、受支持的标准设定（standard setting）。输出是 subgroup-fixed 子空间的笛卡尔正交归一基，可混合母相不可约表示（mix parent irreps）；不采用官网逐 IR 的列约定（per-irrep column convention）。检查依据与缺口见报告。</p></section>
<div class="grid mode-grid">
<section class="panel"><h2>母相 CIF → 本地模式</h2>
<label for="mode-source">母相 CIF 内容</label><textarea id="mode-source" class="mode-input" spellcheck="false" placeholder="粘贴一个独立的普通三维母相 CIF，包含显式对称操作。"></textarea>
<label for="mode-file">或上传 CIF 文件</label><input id="mode-file" type="file" accept=".cif,text/plain">
<small id="mode-input-source">上传保留原始字节；编辑预览后改为 UTF-8 文本。</small>
<div class="row"><div><label for="mode-irrep">IR 标签</label><input id="mode-irrep" type="text" placeholder="例如 GM4-"></div>
<div><label for="mode-opd">OPD 标签</label><input id="mode-opd" type="text" placeholder="例如 P1"></div></div>
<label for="mode-k-values">ISO 独立 k 参数（ISO independent k parameters，按需填写）</label><input id="mode-k-values" type="text" placeholder="空格分隔，例如 1/3">
<small>填写所选 IR 的独立参数；程序不猜测完整 k 分量。CIF 需显式对称操作、空间群编号和完整占位；导出胞最多 256 个原子。</small>
<details><summary>可选振幅：按导出基的列顺序填写</summary>
<label for="mode-amplitudes">位移振幅（Å）</label><input id="mode-amplitudes" type="text" placeholder="有限数值，空格分隔；留空表示全零">
<label for="mode-strain-amplitudes">应变振幅（无量纲）</label><input id="mode-strain-amplitudes" type="text" placeholder="有限数值，空格分隔；留空表示全零">
<p class="hint">先查看 modes.json 的列顺序。振幅个数须与基维数一致；刚体平移数单独列出。</p></details>
<div class="row"><div class="timeout-field"><label for="mode-timeout">每次 ISO 调用超时（秒）</label><input id="mode-timeout" type="number" min="0.1" max="300" step="0.1" value="60"></div><small>最多五次官方调用，总工作流时间可能超过此值。</small></div>
<div class="actions"><button id="mode-run" type="button">计算本地模式</button><small>CIF 上限：1 MiB</small></div>
</section>
<section class="panel"><h2>计算结果与导出</h2>
<div id="mode-status" class="status" role="status">尚未运行模式工作流。</div>
<pre id="mode-summary">此处显示母相、子群、原子数、位移模式数、刚体平移数和应变模式数。</pre>
<div id="mode-downloads" class="downloads"></div>
<details><summary>工作流记录</summary><pre id="mode-record">暂无工作流记录。</pre></details>
<p class="hint">本地子空间检查完成不代表已复现全部 ISODISTORT 功能。复用振幅前，请核对 report.json 与基底约定。</p>
</section></div>
<h2 class="section-heading">Raw official-tool console</h2>
<div class="grid">
<section class="panel"><h2>1. Prepare an official-tool input</h2>
<label for="program">Program</label><select id="program"></select><small id="program-description"></small>
<details class="query-builder" open><summary>Prepare a restricted ISO query</summary>
<div class="row"><div><label for="parent">Parent space group</label><input id="parent" type="number" min="1" max="230" placeholder="1–230"></div>
<div><label for="irrep">IR label</label><input id="irrep" type="text" placeholder="e.g. GM4-"></div></div>
<div class="row"><div><label for="opd">OPD label</label><input id="opd" type="text" placeholder="e.g. P1"></div>
<div><label for="wyckoff">Wyckoff (optional)</label><input id="wyckoff" type="text" maxlength="1" placeholder="e.g. a"></div></div>
<button id="prepare" type="button">Prepare ISO commands</button>
<p class="hint">Shows raw commands for review. Uses your labels as entered; no k parameters or atomic coordinates are inferred. Optional Wyckoff output is an official projection, not a complete normalized CIF mode export.</p>
</details>
<label for="source">Raw input</label><textarea id="source" spellcheck="false" placeholder="Paste trusted input in the selected program's documented format."></textarea>
<div class="row"><div><label for="file">Or upload an input file</label><input type="file" id="file"></div><div class="timeout-field"><label for="timeout">Timeout (s)</label><input type="number" id="timeout" min="0.1" max="300" step="0.1" value="30"></div></div>
<small id="input-source">Text is sent as UTF-8. No coordinates or tolerances are changed.</small>
<div class="actions"><button id="run" disabled>Run official program</button><small>Maximum input: 1 MiB</small></div>
<p class="hint">Use the bundled program manuals for input formats and conventions. Only run inputs you trust: the official programs are not placed in an operating-system sandbox.</p>
</section>
<section class="panel"><h2>2. Inspect the raw result</h2>
<div id="status" class="status" role="status">No computation has been run.</div>
<div class="output-label">Standard output</div><pre id="stdout">Raw program output will appear here.</pre>
<div class="output-label">Standard error</div><pre id="stderr">No output.</pre>
<details><summary>Run evidence and file fingerprints</summary><pre id="evidence">No run record yet.</pre></details>
</section></div>
<footer>Each run has its own working directory and writable copy of the official database files. Inputs, raw output, resource hashes, commands and exit status are saved under <code>runs/&lt;run-id&gt;/</code>.<br>
A zero exit code records process completion only. Scientific input acceptance, numerical accuracy and equivalence to the website require separate validation.</footer>
</main>
<script nonce="__TOKEN__">
const token="__TOKEN__";
const el=id=>document.getElementById(id);
let programs={}, uploaded=null, inputEpoch=0, loadingInput=false, running=false;
let modeUploaded=null, modeEpoch=0, loadingModeInput=false;
function updateRunButton(){el("run").disabled=loadingInput||running||Object.keys(programs).length===0;el("mode-run").disabled=loadingModeInput||running;}
async function call(path, options={}){
 const response=await fetch(path,{...options,headers:{"X-CRIS-Token":token,"Content-Type":"application/json",...(options.headers||{})}});
 const value=await response.json();
 if(!response.ok) throw new Error(value.error||"Request failed");
 return value;
}
function describe(){el("program-description").textContent=programs[el("program").value]?.description||"";}
el("program").addEventListener("change",describe);
el("prepare").addEventListener("click",async()=>{
 const epoch=++inputEpoch;loadingInput=true;updateRunButton();
 try{
  const query=await call("/api/iso-query",{method:"POST",body:JSON.stringify({
   parent:Number(el("parent").value),irrep:el("irrep").value,opd:el("opd").value,wyckoff:el("wyckoff").value||null
  })});
  if(epoch!==inputEpoch)return;
  el("program").value="iso";describe();el("source").value=query.input;uploaded=null;el("file").value="";
  el("input-source").textContent="Prepared official ISO commands. Review the raw input before running.";
  el("status").className="status";el("status").textContent=query.scope;
 }catch(error){if(epoch===inputEpoch){el("status").className="status error";el("status").textContent=error.message;}}
 finally{if(epoch===inputEpoch){loadingInput=false;updateRunButton();}}
});
el("source").addEventListener("input",()=>{++inputEpoch;loadingInput=false;uploaded=null;el("file").value="";el("input-source").textContent="Edited text will be sent as UTF-8.";updateRunButton();});
el("file").addEventListener("change",async()=>{
 const epoch=++inputEpoch;const file=el("file").files[0];uploaded=null;el("source").value="";
 if(!file){loadingInput=false;updateRunButton();return;}
 loadingInput=true;updateRunButton();el("input-source").textContent="Reading "+file.name+"…";
 try{
  if(file.size>1048576)throw new Error("File exceeds 1 MiB.");
  const bytes=new Uint8Array(await file.arrayBuffer());
  if(epoch!==inputEpoch)return;
  let binary="";for(let i=0;i<bytes.length;i++)binary+=String.fromCharCode(bytes[i]);
  uploaded=btoa(binary);el("source").value=new TextDecoder("utf-8").decode(bytes);
  el("input-source").textContent=file.name+" · original file bytes will be sent. Editing the preview switches to UTF-8 text.";
 }catch(error){if(epoch===inputEpoch){el("file").value="";el("input-source").textContent=error.message;}}
 finally{if(epoch===inputEpoch){loadingInput=false;updateRunButton();}}
});
el("run").addEventListener("click",async()=>{
 if(loadingInput||running)return;
 running=true;updateRunButton();el("status").className="status busy";el("status").textContent="Running the selected official program…";
 el("stdout").textContent="";el("stderr").textContent="";el("evidence").textContent="";
 const request={program:el("program").value,timeout:Number(el("timeout").value)};
 if(uploaded!==null)request.input_b64=uploaded;else request.input=el("source").value;
 try{
  const result=await call("/api/run",{method:"POST",body:JSON.stringify(request)});
  const evidence=result.evidence;
  el("status").className="status"+(evidence.status==="completed_uninterpreted"?"":" error");
  el("status").textContent=(evidence.status==="completed_uninterpreted"?"Process completed; scientific result not yet validated.":"Execution status: "+evidence.status)+
   "\nRun: "+evidence.run_id+"\nSaved: "+result.run_directory+(evidence.error?"\n"+evidence.error:"");
  el("stdout").textContent=(result.stdout||"(empty)")+(result.stdout_preview_truncated?"\n[Preview truncated; see stdout.txt for full output.]":"");
  el("stderr").textContent=(result.stderr||"(empty)")+(result.stderr_preview_truncated?"\n[Preview truncated; see stderr.txt for full output.]":"");
  el("evidence").textContent=JSON.stringify(evidence,null,2);
 }catch(error){el("status").className="status error";el("status").textContent=error.message;}
 finally{running=false;updateRunButton();}
});
el("mode-source").addEventListener("input",()=>{
 ++modeEpoch;loadingModeInput=false;modeUploaded=null;el("mode-file").value="";
 el("mode-input-source").textContent="已编辑 CIF 将按 UTF-8 文本提交。";updateRunButton();
});
el("mode-file").addEventListener("change",async()=>{
 const epoch=++modeEpoch;const file=el("mode-file").files[0];modeUploaded=null;el("mode-source").value="";
 if(!file){loadingModeInput=false;updateRunButton();return;}
 loadingModeInput=true;updateRunButton();el("mode-input-source").textContent="正在读取 "+file.name+"…";
 try{
  if(file.size>1048576)throw new Error("CIF 超过 1 MiB。");
  const bytes=new Uint8Array(await file.arrayBuffer());if(epoch!==modeEpoch)return;
  let binary="";for(let i=0;i<bytes.length;i++)binary+=String.fromCharCode(bytes[i]);
  modeUploaded=btoa(binary);el("mode-source").value=new TextDecoder("utf-8").decode(bytes);
  el("mode-input-source").textContent=file.name+" · 将提交原始文件字节。";
 }catch(error){if(epoch===modeEpoch){el("mode-file").value="";el("mode-input-source").textContent=error.message;}}
 finally{if(epoch===modeEpoch){loadingModeInput=false;updateRunButton();}}
});
function modeNumbers(id){
 const value=el(id).value.trim();if(!value)return null;
 const numbers=value.split(/\s+/).map(Number);
 if(numbers.some(value=>!Number.isFinite(value)))throw new Error("振幅须为有限数值（finite numbers），用空格分隔。");
 return numbers;
}
async function downloadMode(workflowId,name){
 const response=await fetch("/api/modes/"+encodeURIComponent(workflowId)+"/"+encodeURIComponent(name),{headers:{"X-CRIS-Token":token}});
 if(!response.ok){const error=await response.json();throw new Error(error.error||"Download failed.");}
 const url=URL.createObjectURL(await response.blob());
 const link=document.createElement("a");link.href=url;link.download=name;document.body.appendChild(link);link.click();link.remove();
 setTimeout(()=>URL.revokeObjectURL(url),1000);
}
el("mode-run").addEventListener("click",async()=>{
 if(loadingModeInput||running)return;
 running=true;updateRunButton();el("mode-status").className="status busy";
 el("mode-status").textContent="正在计算受支持的子群不变子空间…";
 el("mode-summary").textContent="";el("mode-record").textContent="";el("mode-downloads").replaceChildren();
 try{
  const request={irrep:el("mode-irrep").value,opd:el("mode-opd").value,
   k_values:el("mode-k-values").value.trim()?el("mode-k-values").value.trim().split(/\s+/):[],
   amplitudes:modeNumbers("mode-amplitudes"),strain_amplitudes:modeNumbers("mode-strain-amplitudes"),
   timeout:Number(el("mode-timeout").value)};
  if(modeUploaded!==null)request.cif_b64=modeUploaded;else request.cif=el("mode-source").value;
  const result=await call("/api/modes",{method:"POST",body:JSON.stringify(request)});
  const complete=result.status==="validated_local_subspace";
  el("mode-status").className="status"+(complete?"":" error");
  el("mode-status").textContent=(complete?"本地子群不变子空间检查已完成。":"工作流状态："+result.status)+
   "\nIR / OPD: "+request.irrep+" / "+request.opd+"\n工作流："+result.workflow_id+"\n保存位置："+result.directory+
   (result.error?"\n"+(typeof result.error==="string"?result.error:result.error.message||JSON.stringify(result.error)):"");
  const summary=result.summary||{};
  el("mode-summary").textContent=["母相空间群："+(summary.parent??"—"),"子群："+(summary.subgroup??"—")+" "+(summary.subgroup_symbol||""),
   "原子数："+(summary.atoms??"—"),"位移基列数："+(summary.displacement_modes??"—"),"其中刚体平移数："+(summary.rigid_translations??"—"),
   "应变基列数："+(summary.strain_modes??"—"),"基底约定："+(summary.basis_convention||"—")].join("\n");
  el("mode-record").textContent=JSON.stringify(result,null,2);
  if(complete){
   for(const name of ["modes.json","modes.csv","reference.cif","distorted.cif","report.json"]){
    if(!Object.prototype.hasOwnProperty.call(result.files||{},name))continue;
    const button=document.createElement("button");button.type="button";button.textContent="下载 "+name;
    button.addEventListener("click",async()=>{try{await downloadMode(result.workflow_id,name);}catch(error){el("mode-status").className="status error";el("mode-status").textContent=error.message;}});
    el("mode-downloads").appendChild(button);
   }
  }
 }catch(error){el("mode-status").className="status error";el("mode-status").textContent=error.message;}
 finally{running=false;updateRunButton();}
});
call("/api/programs").then(value=>{
 programs=value;for(const name of Object.keys(value)){const option=document.createElement("option");option.value=name;option.textContent=name.toUpperCase();el("program").appendChild(option);}describe();updateRunButton();
}).catch(error=>{el("status").textContent=error.message;el("run").disabled=true;});
</script></body></html>"""
