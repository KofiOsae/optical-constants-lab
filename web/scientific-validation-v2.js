/* v2 Scientific validation: compact, diagnostic-only summary of actual fit data.
   Data absent from results are marked unavailable, never inferred. */
(()=>{'use strict';function init(){
const host=document.getElementById('diagnostics')||document.querySelector('.diagnostics')||document.getElementById('results')||document.body;if(!host)return;
const section=document.createElement('details');section.id='v2-science-validation';section.style.cssText='margin:8px 0;border:1px solid #4a647c;border-radius:8px;padding:10px;';
section.innerHTML='<summary style="font-weight:600;cursor:pointer">Scientific validation · v2</summary><p style="font-size:12px;opacity:.85">Compact summary of diagnostics available from the current completed fit. No additional fitting is launched.</p><div id="v2-validation-tabs" style="display:flex;flex-wrap:wrap;gap:5px"></div><div id="v2-validation-content" style="white-space:pre-wrap;line-height:1.55;font-size:13px;padding:8px"></div>';
host.append(section);
const tabs=[['fit','Fit quality'],['reliability','Reliability'],['physics','Physics'],['prediction','Prediction'],['assisted','Assisted assessment']];let active='fit';
const get=()=>typeof window.oclValidationResult==='function'?window.oclValidationResult():null;
const fmt=v=>Number.isFinite(Number(v))?Number(v).toPrecision(5):'Not reported';
function assisted(d){
 if(!d)return 'NOT ASSESSED — no completed fit is available. Run a fit and reopen this panel.';
 const sens=d.sensitivity||{},preds=Array.isArray(d.predictions)?d.predictions:[],out=[];
 const add=(status,label,reason)=>out.push(status+' · '+label+': '+reason);
 const success=d.success??d.converged;
 if(success===false)add('FAIL','Optimizer','Solver did not report convergence; inspect status and evaluation budget.');
 else if(success===true)add('PASS (numerical only)','Optimizer','Convergence was reported; this does not establish model adequacy or unique optical constants.');
 else add('NOT ASSESSED','Optimizer','Convergence flag not available in this result.');
 const weighted=Number(d.weighted_rmse);
 if(d.weighted_rmse!==undefined&&Number.isFinite(weighted))add('REVIEW','Spectral fit','Weighted RMSE '+weighted.toPrecision(5)+'. Compare residuals against the independently known R/T noise level and inspect structured residuals; no universal RMSE pass threshold.');
 else add('NOT ASSESSED','Spectral fit','No weighted RMSE returned.');
 const rank=Number(sens.rank),count=Number(sens.parameters);
 if(sens.rank!=null&&sens.parameters!=null&&Number.isFinite(rank)&&Number.isFinite(count)&&count>0){
   if(rank<count)add('FAIL (local identifiability)','Jacobian rank',rank+'/'+count+' free parameter directions; at least '+(count-rank)+' direction(s) not independently constrained at numerical tolerance.');
   else if(rank===count)add('PASS (local rank only)','Jacobian rank',rank+'/'+count+'; full numerical rank does not prove global uniqueness.');
   else add('REVIEW','Jacobian rank','Reported rank exceeds parameter count: check solver diagnostics.');
 }else add('NOT ASSESSED','Jacobian rank','Rank and/or free-parameter count absent; do not infer full rank from rank alone.');
 const kappa=Number(sens.condition_number);
 if(sens.condition_number!=null&&Number.isFinite(kappa)&&kappa>=1){
   add('REVIEW','Conditioning','κ='+kappa.toPrecision(5)+' (reported Jacobian). Interpret only with parameter scaling, weighting, singular values and perturbation tests; no universal pass/fail cutoff.');
 }else add('NOT ASSESSED','Conditioning','Condition number absent or invalid.');
 const starts=d.stages?.start_results||d.multistart?.runs||d.multistart?.start_results||[];
 const attempts=Number(d.multistart?.attempts)||starts.length;
 if(Array.isArray(starts)&&starts.length>1){
   const vectors=starts.every(x=>Array.isArray(x.coefficients)||Array.isArray(x.parameters));
   add('REVIEW','Multistart',starts.length+' recorded starts; '+(vectors?'compare their final n,k or parameter vectors at similarly low costs.':'fitted vectors unavailable, so cost agreement alone cannot prove identical recovered material.'));
 }else if(attempts>1)add('NOT ASSESSED','Multistart',attempts+' starts reported but complete fitted vectors not available to this panel; check the individual-start exports.');
 else add('NOT ASSESSED','Multistart','No independent starting solutions reported.');
 const k=Array.isArray(d.k)?d.k:[],n=Array.isArray(d.n)?d.n:[];
 if(k.length&&n.length){
   const minK=Math.min(...k),bad=n.some(v=>!Number.isFinite(Number(v))||Number(v)<=0)||k.some(v=>!Number.isFinite(Number(v)));
   if(bad)add('FAIL','Physical constants','Non-finite n/k or non-positive n detected; verify the recovered optical model.');
   else if(minK< -1e-8)add('FAIL','Passivity','Negative extinction coefficient: minimum k='+minK.toPrecision(5)+' under n+i k convention.');
   else add('PASS (basic physics)','Passivity','Finite positive n and k≥0; dispersion, causality, materials and stack consistency still need independent checks.');
 }else add('NOT ASSESSED','Passivity','Recovered n/k arrays unavailable.');
 if(preds.length&&preds.every(p=>Array.isArray(p.R_fit)&&Array.isArray(p.T_fit))){
   const excess=preds.flatMap(p=>p.R_fit.map((v,i)=>Number(v)+Number(p.T_fit[i])-1)).filter(Number.isFinite);
   const maxExcess=excess.length?Math.max(...excess):NaN;
   if(Number.isFinite(maxExcess)&&maxExcess>1e-5)add('REVIEW','Energy balance','R+T exceeds 1 in predicted spectra. Confirm normal flux normalization and substrate/collection assumptions.');
   else if(Number.isFinite(maxExcess))add('PASS (conditional)','Energy balance','No R+T>1 violation detected for normalized passive-stack assumptions.');
 }else add('NOT ASSESSED','Energy balance','Predicted R/T arrays not available for check.');
 add('NOT ASSESSED','Material accuracy','Independent same-material n/k or analytical benchmark comparison must be performed; a known-layer Au file is not an SiO₂ reference.');
 add('NOT ASSESSED','Prediction and robustness','Independent withheld data and noise/thickness/window perturbations must be separately run and assessed.');
 const fail=out.filter(x=>x.startsWith('FAIL')).length;
 return 'ASSISTED SCIENTIFIC FIT REVIEW — '+(fail?'CHECK FAILURES':'PROVISIONAL; NO AUTOMATIC SCIENTIFIC CERTIFICATION')+'\n'+out.join('\n')+'\n\nSuggested next actions: verify layer/reference provenance; compare R/T residuals with actual noise; examine fitted n/k against independent truth; export parameter vectors for all starts; perturb wavelength and thickness; predict truly held-out data. See Guided tutorial → Fit reliability for interpretation. PASS labels denote only the stated local numerical check.';
}
function text(d,tab){if(tab==='assisted')return assisted(d);if(!d)return 'No completed fit is currently available for this workflow. Fit history CSV records are summaries of earlier runs, not full parameter/Jacobian data. Run a fresh fit, then reopen this panel.';
if(tab==='fit'){const preds=d.predictions||[];return 'Global weighted RMSE: '+fmt(d.weighted_rmse)+'\n'+preds.map((p,i)=>'Sample '+(i+1)+' R/T RMSE: '+fmt(p.rmse)).join('\n')+'\nFor systematic wavelength-dependent residuals, inspect the existing Plotly residual chart. A low scalar RMSE alone does not establish correct optical constants.'}
if(tab==='reliability'){const sens=d.sensitivity||{};const lines=['Jacobian rank: '+(sens.rank??'Not reported'),'Fitted parameters: '+(sens.parameters??'Not reported'),'Jacobian condition number: '+fmt(sens.condition_number),'Local standard errors: '+(Array.isArray(sens.parameter_standard_errors_local)?sens.parameter_standard_errors_local.map(fmt).join(', '):'Not reported')];if(sens.identifiability_warnings?.length)lines.push('Warnings: '+sens.identifiability_warnings.join(' | '));const starts=d.stages?.start_results||[];const m=d.multistart||{};const ambiguity=d.ambiguity||{};if(starts.length){const scores=starts.map(x=>Number(x.weighted_rmse)).filter(Number.isFinite);lines.push('Multistart runs: '+starts.length+'; weighted RMSE range: '+(scores.length?Math.min(...scores).toPrecision(4)+'–'+Math.max(...scores).toPrecision(4):'Not reported'))}else if(m.attempts){lines.push('Multistart attempts: '+m.attempts+'; distinct near-optimal branches: '+(ambiguity.near_optimal_branches??'Not reported')+'; interpretation: '+(ambiguity.status||'Unreported'))}else lines.push('Multi-start stability: not measured in this run.');lines.push('Thickness-profile and wavelength-window stability need separate refits; one fitted solution cannot provide those tests.');return lines.join('\n')}
if(tab==='physics'){const w=d.wavelength_nm||[],n=d.n||[],k=d.k||[];let minK=k.length?Math.min(...k):null;const preds=d.predictions||[];let maxExcess=null;for(const p of preds){const R=p.R_fit,T=p.T_fit;if(Array.isArray(R)&&Array.isArray(T)){for(let j=0;j<Math.min(R.length,T.length);j++){const excess=Number(R[j])+Number(T[j])-1;if(Number.isFinite(excess))maxExcess=Math.max(maxExcess??-Infinity,excess)}}}return 'Minimum recovered extinction k: '+(minK===null?'Not reported':fmt(minK))+'\nMaximum (R_fit + T_fit − 1), if available: '+(maxExcess===null?'Not reported':fmt(maxExcess))+'\nNote: these checks assume normalized flux and passive stack conventions. KK/SSKK finite-band comparisons are separate diagnostics and are not proofs of causality. Absorbing substrates need special interpretation.'}
return 'Out-of-sample predictive validation must be based on a deliberately excluded measurement, or a separately generated reference forward solution.\n'+(d.holdout?'Holdout information present in fit result. Inspect its explicit prediction metrics.':'No explicit held-out validation confirmed from this result.')+'\nModel comparison (AICc/BIC), spectral-window stability, and multi-start consistency are not automatically calculated by this summary.'}
function render(){const d=get();const area=section.querySelector('#v2-validation-content');if(!area)return;try{
 const report=text(d,active);
 if(active==='assisted'){
  area.replaceChildren();
  area.style.whiteSpace='normal';
  const lines=report.split('\n');
  const heading=document.createElement('p');heading.textContent=lines.shift()||'Assisted assessment';heading.style.cssText='font-weight:700;margin:0 0 12px;color:#e4f0ff';area.append(heading);
  const list=document.createElement('div');list.style.cssText='display:grid;gap:8px';
  for(const line of lines){
   if(!line.trim())continue;
   const sep=line.indexOf(' · '),colon=line.indexOf(': ',sep+3);
   if(sep<0||colon<0){const p=document.createElement('p');p.textContent=line;p.style.cssText='margin:8px 0;line-height:1.55';area.append(p);continue}
   const status=line.slice(0,sep),label=line.slice(sep+3,colon),reason=line.slice(colon+2);
   const item=document.createElement('div');item.style.cssText='border:1px solid #405a72;border-radius:7px;padding:10px;background:rgba(12,28,48,.46)';
   const h=document.createElement('div');h.style.cssText='font-weight:650;margin-bottom:5px';
   const tag=document.createElement('span');tag.textContent=status;tag.style.cssText='display:inline-block;font-size:11px;font-weight:700;letter-spacing:.02em;padding:2px 5px;margin-right:8px;border:1px solid #718ba6;border-radius:4px';
   const name=document.createElement('span');name.textContent=label;h.append(tag,name);
   const detail=document.createElement('div');detail.textContent=reason;detail.style.cssText='font-size:13px;line-height:1.55;overflow-wrap:anywhere';
   item.append(h,detail);list.append(item);
  }
  area.insertBefore(list,heading.nextSibling);
 }else{area.style.whiteSpace='pre-wrap';area.textContent=report;}
}catch(error){area.textContent='Validation summary unavailable: '+error.message+'; see browser console for details.';}section.querySelectorAll('button').forEach(b=>{b.style.background=b.dataset.tab===active?'#2879a4':'#233b56'})}
for(const [key,label] of tabs){const b=document.createElement('button');b.type='button';b.dataset.tab=key;b.textContent=label;b.style.cssText='padding:6px 10px;border:1px solid #6489a6;border-radius:6px;color:#fff';b.onclick=()=>{active=key;render()};section.querySelector('#v2-validation-tabs').append(b)}
section.addEventListener('toggle',()=>{if(section.open)render()});document.addEventListener('ocl-fit-complete',render);render();
}if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',init);else init()})();
/* Cross-workflow collapsible groups. Move existing DOM controls without cloning or
   changing element IDs: event listeners, file inputs, results and save/load state survive. */
(()=>{'use strict';
function init(){
 const one=id=>document.getElementById(id);
 function group(parent,items,label,open=true){
  const nodes=items.filter(Boolean).filter(n=>n.parentElement===parent);
  if(!nodes.length)return null;
  const details=document.createElement('details');details.className='v2-collapse-group';details.open=open;
  details.style.cssText='border:1px solid #486078;border-radius:8px;padding:9px;margin:9px 0;min-width:0';
  const summary=document.createElement('summary');summary.textContent=label;
  summary.style.cssText='font-weight:650;cursor:pointer;padding:4px 2px';
  parent.insertBefore(details,nodes[0]);details.append(summary,...nodes);
  details.addEventListener('toggle',()=>{
   if(!details.open||!window.Plotly?.Plots?.resize)return;
   requestAnimationFrame(()=>details.querySelectorAll('.js-plotly-plot').forEach(plot=>{
    try{window.Plotly.Plots.resize(plot)}catch(_){}
   }));
  });
  return details;
 }
 const isMulti=!!one('controls')&&!!one('diagnostics')&&!one('single_upload');
 if(!isMulti){
  const single=one('single_upload'),joint=one('joint_upload');
  if(single&&joint&&single.parentElement===joint.parentElement)
   group(single.parentElement,[single,joint],'Measured spectra · import and preview',true);
  const ref=one('reference_ri_file')?.closest('.panel');
  if(ref&&ref.querySelector('h3')?.textContent.includes('Reference optical constants comparison')){
   const title=ref.querySelector('h3');title.remove();
   const d=document.createElement('details');d.open=true;d.className='v2-collapse-group';
   d.style.cssText='min-width:0';
   const summary=document.createElement('summary');summary.textContent='Reference optical constants comparison';summary.style.cssText='cursor:pointer;font-weight:700';
   while(ref.firstChild)d.append(ref.firstChild);
   d.prepend(summary);ref.append(d);
  }
  const diagnostics=document.querySelector('aside.diagnostics');
  if(diagnostics){
   group(diagnostics,[one('metrics'),one('diagnostic_notes')],'Fit metrics and diagnostic interpretation',true);
   const history=one('history_view')?.closest('section');
   if(history&&history.parentElement===diagnostics){
    const heading=history.querySelector('h3');if(heading)heading.remove();
    group(diagnostics,[history],'Fit history · RMSE and starting solutions',true);
   }
   const review=one('review');
   if(review&&review.parentElement===diagnostics)group(diagnostics,[review],'Phase I inspection',true);
  }
 }else{
  const controls=one('controls');
  if(controls){
   for(const fieldset of controls.querySelectorAll(':scope > fieldset')){
    if(fieldset.closest('details'))continue;
    const label=fieldset.querySelector(':scope > legend')?.textContent?.trim()||'Data & settings';
    group(controls,[fieldset],label,true);
   }
  }
  const diagnostics=one('diagnostics');
  if(diagnostics){
   group(diagnostics,[one('metrics'),one('rmsePlot')],'Fit metrics and RMSE',true);
   group(diagnostics,[one('multistartPlot'),one('branchsummary')],'Multistart and solution reliability',true);
  }
 }
}
if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',init);else init();
})();
