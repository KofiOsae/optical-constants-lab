/* Private, opt-in quality register for V2. Stores only local counters and categorized failures.
   No telemetry is transmitted and no raw spectra, filenames, addresses, identities,
   fitting parameters, or measurements are recorded. */
(()=>{'use strict';
function start(){
 const nav=document.querySelector('nav.main-tabs');if(!nav||document.getElementById('ocl-quality'))return;
 const KEY='ocl-v2-local-quality-v1';
 const read=()=>{try{const v=JSON.parse(localStorage.getItem(KEY)||'{}');return {enabled:!!v.enabled,counts:v.counts||{},failures:v.failures||{}}}catch{return {enabled:false,counts:{},failures:{}}}};
 let state=read();
 const save=()=>{try{localStorage.setItem(KEY,JSON.stringify(state))}catch{}};
 function count(key){if(!state.enabled)return;state.counts[key]=(state.counts[key]||0)+1;save();render()}
 const button=document.createElement('button');button.id='ocl-quality';button.type='button';button.textContent='Quality & usage';nav.append(button);
 const modal=document.createElement('dialog');modal.id='ocl-quality-dialog';modal.style.cssText='width:min(540px,94vw);max-height:78vh;overflow:auto;background:#14253a;color:#eef5fb;border:1px solid #6186a5;border-radius:12px;padding:18px';
 modal.innerHTML='<h2>Local laboratory quality register</h2><p>Optional, device-local only. No data is transmitted to the laboratory. This is not a count of all users. Raw data, filenames, fit parameters and personal identifiers are never collected automatically.</p><label><input id="ocl-quality-enable" type="checkbox"> Record anonymized interaction counts on this device</label><div id="ocl-quality-counts"></div><p><strong>Log a failure category</strong> (no free text or personal data)</p><select id="ocl-quality-category"><option value="fit">Fit/solver</option><option value="plot">Plots/visualization</option><option value="export">Exports</option><option value="project">Save/Load</option><option value="diagnostic">Diagnostics</option><option value="controls">Controls/layout</option></select><button id="ocl-quality-log" type="button">Record category</button><div id="ocl-quality-failures"></div><div style="display:flex;flex-wrap:wrap;gap:8px"><button id="ocl-quality-export" type="button">Export aggregate JSON</button><button id="ocl-quality-clear" type="button">Clear local records</button><button id="ocl-quality-close" type="button">Close</button></div>';
 document.body.append(modal);
 const $=id=>modal.querySelector('#'+id);
 function render(){
  $('ocl-quality-enable').checked=state.enabled;
  $('ocl-quality-counts').textContent='Interaction counts: '+(Object.entries(state.counts).map(([k,v])=>k+': '+v).join(' · ')||'none');
  $('ocl-quality-failures').textContent='Logged failure categories: '+(Object.entries(state.failures).map(([k,v])=>k+': '+v).join(' · ')||'none');
 }
 button.onclick=()=>{render();modal.showModal()};
 $('ocl-quality-close').onclick=()=>modal.close();
 $('ocl-quality-enable').onchange=e=>{state.enabled=e.target.checked;save();render()};
 $('ocl-quality-log').onclick=()=>{if(!state.enabled){alert('Enable local recording first.');return}const k=$('ocl-quality-category').value;state.failures[k]=(state.failures[k]||0)+1;save();render()};
 $('ocl-quality-clear').onclick=()=>{if(!confirm('Clear all local quality statistics?'))return;state={enabled:false,counts:{},failures:{}};localStorage.removeItem(KEY);render()};
 $('ocl-quality-export').onclick=()=>{const snapshot={schema:'ocl-v2-local-quality-v1',counts:state.counts,failures:state.failures,scope:'one browser only; no network telemetry'};const url=URL.createObjectURL(new Blob([JSON.stringify(snapshot,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='optical_lab_anonymous_quality_counts.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)};
 document.addEventListener('click',e=>{const id=e.target.closest('button')?.id;if(id==='run')count('single_joint_fit_requests');if(id==='download')count('optical_csv_requests');if(id==='download_model_xlsx')count('model_workbook_requests');if(id==='saveWorkspaceProject')count('project_save_requests')},true);
 const frame=document.getElementById('multilayer_frame');
 frame?.addEventListener('load',()=>{try{frame.contentDocument?.addEventListener('click',e=>{const id=e.target.closest('button')?.id;if(id==='run')count('multilayer_fit_requests');if(id==='exportModelExcel')count('multilayer_workbook_requests')},true)}catch{}});
 render();
}
if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',start);else start();
})();
