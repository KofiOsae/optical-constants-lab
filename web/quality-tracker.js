/* Anonymous aggregate usage dashboard. No cookies, IDs, fingerprints, IP addresses,
   filenames, fitted coefficients or uploaded spectra are transmitted by this script.
   "Visits" count page loads, not distinct visitors. */
(()=>{'use strict';
function start(){
 const slot=document.getElementById('ocl-stats-slot')||document.querySelector('nav.main-tabs');
 if(!slot||document.getElementById('ocl-quality'))return;
 const button=document.createElement('button');
 button.id='ocl-quality';button.type='button';button.textContent='User statistics';
 button.title='Anonymous daily visit and fitting activity (not unique users)';
 slot.append(button);
 const dialog=document.createElement('dialog');
 dialog.id='ocl-quality-dialog';
 dialog.style.cssText='width:min(760px,96vw);max-height:85vh;overflow:auto;background:#14253a;color:#eef5fb;border:1px solid #6186a5;border-radius:12px;padding:18px';
 dialog.innerHTML='<h2>User statistics · anonymous aggregate activity</h2><p style="line-height:1.5">Counts page visits, completed fitting requests and exports over time. No visitor identifiers, cookies, IP addresses or uploaded spectra are recorded in the usage database. <strong>Visits are not unique users.</strong> Automated traffic can affect counts.</p><label for="ocl-usage-range">Period</label><select id="ocl-usage-range"><option value="7">Past 7 days</option><option value="30" selected>Past 30 days</option><option value="365">Past 12 months</option></select><div id="ocl-usage-status" role="status" aria-live="polite"></div><div id="ocl-usage-totals" style="display:flex;flex-wrap:wrap;gap:12px;margin:10px 0"></div><div id="ocl-usage-chart" style="height:290px"></div><div id="ocl-usage-workflows"></div><p style="font-size:12px;opacity:.8">A successful fitting API response is counted once. CSV exports and project saves count browser actions. Totals are approximate; no person-level frequency or unique-user tracking occurs.</p><button id="ocl-quality-close" type="button">Close</button>';
 document.body.append(dialog);
 const $=id=>dialog.querySelector('#'+id);
 function send(event){
   try{fetch('/api/usage/event',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({event}),credentials:'omit',cache:'no-store',keepalive:true,referrerPolicy:'no-referrer'}).catch(()=>{})}catch(_){}
 }
 let visitSent=false; if(!visitSent){visitSent=true;send('visit')}
 async function render(){
  $('ocl-usage-status').textContent='Loading aggregate usage…';
  try{
   const response=await fetch('/api/usage/summary?days='+$('ocl-usage-range').value,{cache:'no-store',credentials:'omit'});
   if(!response.ok)throw Error('Statistics unavailable (HTTP '+response.status+')');
   const data=await response.json();
   if(!data.enabled){$('ocl-usage-status').textContent=data.reason||'Persistent aggregate storage is not connected.';$('ocl-usage-totals').replaceChildren();$('ocl-usage-chart').replaceChildren();$('ocl-usage-workflows').textContent='No data have been stored.';return}
   $('ocl-usage-status').textContent='UTC daily aggregates · '+data.days+' days';
   const totals=data.totals||{};
   const cards=[['Visits (page loads)',totals.visit||0],['Completed fit requests',(totals.fit_single||0)+(totals.fit_joint||0)+(totals.fit_multilayer||0)],['Fit requests per visit',totals.visit?(((totals.fit_single||0)+(totals.fit_joint||0)+(totals.fit_multilayer||0))/totals.visit).toFixed(2):'—'],['CSV exports',totals.export_csv||0],['Excel exports',totals.export_excel||0]];
   $('ocl-usage-totals').replaceChildren(...cards.map(([title,value])=>{const el=document.createElement('div');el.style.cssText='padding:10px;border:1px solid #47627e;border-radius:7px;min-width:126px';const name=document.createElement('div');name.textContent=title;name.style.cssText='font-size:11px;opacity:.8';const number=document.createElement('strong');number.textContent=String(value);number.style.fontSize='22px';el.append(name,number);return el}));
   const labels=data.daily.map(p=>p.date),visits=data.daily.map(p=>p.visit||0),fits=data.daily.map(p=>(p.fit_single||0)+(p.fit_joint||0)+(p.fit_multilayer||0));
   if(window.Plotly){window.Plotly.react($('ocl-usage-chart'),[{x:labels,y:visits,name:'Visits',type:'scatter',mode:'lines+markers'},{x:labels,y:fits,name:'Fitted requests',type:'scatter',mode:'lines+markers'}],{title:{text:'Daily activity (UTC)',font:{size:13}},paper_bgcolor:'#14253a',plot_bgcolor:'#14253a',font:{color:'#e7effa',size:11},margin:{l:40,r:12,t:38,b:47},legend:{orientation:'h',y:-.22},xaxis:{automargin:true},yaxis:{rangemode:'tozero',title:'Count',automargin:true}},{responsive:true,displayModeBar:false})}
   $('ocl-usage-workflows').textContent='Completed fitting requests — Single: '+(totals.fit_single||0)+' · Joint: '+(totals.fit_joint||0)+' · Multilayer: '+(totals.fit_multilayer||0)+' · Project saves: '+(totals.project_save||0);
  }catch(e){$('ocl-usage-status').textContent='Statistics unavailable: '+e.message}
 }
 button.onclick=()=>{dialog.showModal();render()};
 $('ocl-quality-close').onclick=()=>dialog.close();
 $('ocl-usage-range').onchange=render;
 document.addEventListener('click',e=>{
  const id=e.target.closest('button')?.id;
  if(id==='download'||id==='export')send('export_csv');
  if(id==='ocl-top-save')send('project_save');
 },true);
 const frame=document.getElementById('multilayer_frame');
 frame?.addEventListener('load',()=>{
  try{frame.contentDocument?.addEventListener('click',e=>{const id=e.target.closest('button')?.id;if(id==='export')send('export_csv')},true)}catch(_){}
 });
}
if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',start);else start();
})();
