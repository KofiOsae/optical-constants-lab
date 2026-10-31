/* Optional Plotly export helper for v2. Existing Plotly plots stay fully interactive. */
(function(){
'use strict';
function init(){
 if(typeof Plotly==='undefined')return;
 const style=document.createElement('style');
 style.textContent='#ocl-saveplot{position:fixed;bottom:16px;right:16px;z-index:10000;background:#235d82;color:white;border:1px solid #7bb6d6;border-radius:10px;padding:10px 14px;font-weight:600}#ocl-saveplot-dialog{width:min(560px,94vw);max-height:85vh;overflow:auto;background:#15243a;color:#edf5ff;border:1px solid #668daa;border-radius:12px;padding:20px}#ocl-saveplot-dialog::backdrop{background:#050b17b3}#ocl-saveplot-dialog button,#ocl-saveplot-dialog select{background:#244869;color:white;padding:9px;margin:5px;border:1px solid #8aaec5;border-radius:7px}';
 document.head.append(style);
 const open=document.createElement('button');open.type='button';open.id='ocl-saveplot';open.textContent='Save a plot';
 const dialog=document.createElement('dialog');dialog.id='ocl-saveplot-dialog';
 dialog.innerHTML='<h2>Save Plotly figure</h2><p>Choose any visible scientific plot. The existing interactive display, hover values, axes and numerical solver remain unchanged.</p><label>Plot <select id="ocl-plot-select"></select></label><p><button type="button" id="ocl-save-svg">Save SVG (vector)</button><button type="button" id="ocl-save-png">Save PNG (3000 px wide)</button><button type="button" id="ocl-save-csv">Save plotted data CSV</button><button type="button" id="ocl-save-close">Close</button></p><p id="ocl-save-status" aria-live="polite"></p>';
 document.body.append(open,dialog);
 const get=id=>document.getElementById(id);
 const status=msg=>get('ocl-save-status').textContent=msg;
 const plots=()=>[...document.querySelectorAll('.js-plotly-plot')].filter(p=>p.id&&Array.isArray(p.data)&&p.data.length&&p.getClientRects().length&&p.offsetWidth>0);
 function choose(){const all=plots();const select=get('ocl-plot-select');select.replaceChildren();all.forEach(p=>{const label=(p.layout&&p.layout.title&&p.layout.title.text)||p.id;select.add(new Option(String(label).replace(/<[^>]*>/g,''),p.id))});status(all.length?all.length+' plots available for export.':'No plotted data yet. Preview or fit your spectra first.')}
 function plot(){const id=get('ocl-plot-select').value;const p=get(id);if(!p||!p.data)throw Error('Select a visible scientific plot');return p}
 function download(url,name){const a=document.createElement('a');a.href=url;a.download=name;document.body.append(a);a.click();a.remove()}
 function base(){return 'optical_lab_'+plot().id.replace(/[^a-z0-9_-]/gi,'_')}
 async function run(task){try{status('Preparing export…');await task();status('Export prepared; check browser downloads.')}catch(err){status('Export failed: '+err.message)}}
 open.onclick=()=>{choose();dialog.showModal()};
 get('ocl-save-close').onclick=()=>dialog.close();
 get('ocl-save-svg').onclick=()=>run(async()=>{const p=plot();const url=await Plotly.toImage(p,{format:'svg',width:1100,height:700});download(url,base()+'.svg')});
 get('ocl-save-png').onclick=()=>run(async()=>{const p=plot();const url=await Plotly.toImage(p,{format:'png',width:1500,height:950,scale:2});download(url,base()+'.png')});
 get('ocl-save-csv').onclick=()=>run(async()=>{const p=plot();const rows=['trace,x,y'];function esc(x){return '"'+String(x==null?'':x).replaceAll('"','""')+'"'}for(const t of p.data){const xs=t.x||[],ys=t.y||[];for(let i=0;i<Math.min(xs.length,ys.length);i++)rows.push([esc(t.name||'trace'),esc(xs[i]),esc(ys[i])].join(','))}const url=URL.createObjectURL(new Blob([rows.join('\n')],{type:'text/csv'}));download(url,base()+'.csv');setTimeout(()=>URL.revokeObjectURL(url),1000)});
}
if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',init);else init();
})();