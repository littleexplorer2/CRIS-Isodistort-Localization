PAGE = r"""<!doctype html>
<html lang="en">
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
@media(max-width:800px){.grid{grid-template-columns:1fr}header{padding:22px 24px}.local{display:none}}
</style>
</head><body>
<header><div><h1>CRIS REA</h1><p>Official ISOTROPY programs · local execution · traceable runs</p></div><span class="local">● &nbsp; Loopback only</span></header>
<main>
<section class="target"><div class="eyebrow">Research priority</div>
<h2>Parent structure + k / IR / OPD → subgroups and distortion modes</h2>
<p>The current console runs the official tools with your raw input. A verified parent-structure workflow is planned; this console does not yet reproduce the ISODISTORT website methods.</p></section>
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
function updateRunButton(){el("run").disabled=loadingInput||running||Object.keys(programs).length===0;}
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
call("/api/programs").then(value=>{
 programs=value;for(const name of Object.keys(value)){const option=document.createElement("option");option.value=name;option.textContent=name.toUpperCase();el("program").appendChild(option);}describe();updateRunButton();
}).catch(error=>{el("status").textContent=error.message;el("run").disabled=true;});
</script></body></html>"""
