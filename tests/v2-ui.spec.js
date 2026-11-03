const {test,expect}=require('@playwright/test');
const base='http://127.0.0.1:8765';
test('Single V2 input and navigation parity',async({page})=>{
 await page.goto(base);
 await expect(page.locator('#workflow')).toBeVisible();
 await expect(page.locator('#ocl-top-save')).toBeVisible();
 await expect(page.locator('#ocl-top-load')).toBeVisible();
 await expect(page.locator('#ocl-kk-control')).toHaveCount(1);
 await page.locator('#file').setInputFiles('examples/v2/transparent_film_120nm.csv');
 await expect(page.locator('#thickness')).toHaveValue('120');
 await expect(page.locator('#reference_ri_file')).toHaveCount(1);
 await expect(page.locator('#ocl-reference-inputs')).toHaveCount(1);
 await page.locator('#nav_tutorial').click();
 await expect(page.locator('#tutorial_main #ocl-examples')).toHaveCount(1);
});
test('Joint inputs expose inferred thickness, shared reference and progressive controls',async({page})=>{
 await page.goto(base);
 await page.locator('#workflow').selectOption('joint');
 await page.locator('#files').setInputFiles(['examples/v2/transparent_film_120nm.csv','examples/v2/transparent_film_240nm.csv']);
 await expect(page.locator('#specimens .file-card')).toHaveCount(2);
 await expect(page.locator('#ocl-joint-stage')).toHaveCount(1);
 await expect(page.locator('#reference_ri_file')).toHaveCount(1);
});
test('Embedded Multilayer filename thickness and diagnostic controls',async({page})=>{
 await page.goto(base);
 await page.locator('#workflow').selectOption('multilayer');
 const frame=page.frameLocator('#multilayer_frame');
 await expect(frame.locator('#run')).toBeVisible();
 await frame.locator('.file').first().setInputFiles('examples/v2/transparent_film_120nm.csv');
 await expect(frame.locator('.thickness').first()).toHaveValue('120');
 await expect(frame.locator('#startsExport')).toHaveCount(1);
 await expect(frame.locator('#multistartPlot')).toHaveCount(1);
 await expect(frame.locator('#reference')).toHaveCount(1);
});

test('Visual layout regression: no duplicated visible embedded Multilayer guides',async({page})=>{
 await page.goto(base);
 await page.locator('#workflow').selectOption('multilayer');
 const frame=page.frameLocator('#multilayer_frame');
 await expect(frame.locator('#v2StackDesigner')).toBeVisible();
 await expect(frame.locator('#widthControls')).toBeHidden();
 await expect(frame.locator('h1')).toBeHidden();
 await expect(frame.getByText('Guided example datasets',{exact:false})).toBeHidden();
 await expect(frame.getByText('Measurement suitability — read before fitting experimental R/T',{exact:false})).toBeHidden();
 await expect(frame.locator('#stageControls')).toHaveCount(1);
 await expect(frame.locator('#reference')).toHaveCount(1);
});
test('Single and Joint do not show a duplicate measurement suitability shortcut',async({page})=>{
 await page.goto(base);
 await expect(page.getByRole('button',{name:'Measurement suitability · assumptions'})).toBeHidden();
 await page.locator('#workflow').selectOption('joint');
 await expect(page.getByRole('button',{name:'Measurement suitability · assumptions'})).toBeHidden();
});

test('Shared column sliders govern embedded Multilayer and survive workflow switches',async({page})=>{
 await page.goto(base);
 await page.locator('#workflow').selectOption('multilayer');
 const frame=page.frameLocator('#multilayer_frame');
 await expect(frame.locator('#layout')).toBeVisible();
 await page.locator('#singleColumnWidth').fill('34');
 await page.locator('#diagnosticColumnWidth').fill('20');
 await expect(frame.locator('html')).toHaveCSS('--settings-width','34%');
 await expect(frame.locator('html')).toHaveCSS('--diagnostics-width','20%');
 await page.locator('#workflow').selectOption('single');
 await page.locator('#workflow').selectOption('multilayer');
 await expect(frame.locator('html')).toHaveCSS('--settings-width','34%');
 await page.locator('#resetColumnWidths').click();
 await expect(frame.locator('html')).toHaveCSS('--settings-width','27%');
 await expect(frame.locator('html')).toHaveCSS('--diagnostics-width','25%');
});

test('Multilayer controls follow shared input sequence and adapt to the selected model',async({page})=>{
 await page.goto(base);
 await page.locator('#workflow').selectOption('multilayer');
 const frame=page.frameLocator('#multilayer_frame');
 const order=await frame.locator('#controls').evaluate(el=>{
  const ids=['samples','reference','model','stageControls','diagnostic'];
  return ids.map(id=>Array.from(el.querySelectorAll('*')).indexOf(el.querySelector('#'+id)));
 });
 expect(order.every((v,i)=>v>=0&&(i===0||v>order[i-1]))).toBeTruthy();
 await frame.locator('#model').locator('xpath=ancestor::details[1]/summary').click();
 await frame.locator('#model').selectOption('spline');
 await expect(frame.locator('#basis').locator('..')).toBeVisible();
 await expect(frame.locator('#osc').locator('..')).toBeHidden();
 await frame.locator('#model').selectOption('lorentz');
 await expect(frame.locator('#osc').locator('..')).toBeVisible();
 await expect(frame.locator('#basis').locator('..')).toBeHidden();
});

test('Multilayer uses only one R/T graph on preview and retains roughness entry',async({page})=>{
 await page.goto(base);await page.locator('#workflow').selectOption('multilayer');
 const f=page.frameLocator('#multilayer_frame');
 await f.locator('.file').first().setInputFiles('examples/v2/transparent_film_120nm.csv');
 await expect(f.locator('.interfaces').first()).toBeAttached();
 await expect(f.locator('#saveProject')).toBeHidden();
 await expect(f.locator('#loadProject').locator('..')).toBeHidden();
 await expect(f.locator('#preview')).toHaveCount(0);
 await f.locator('#previewButton').click();
 await expect(f.locator('#spectra.js-plotly-plot')).toHaveCount(1,{timeout:15000});
 await expect(f.locator('#viewSample option')).toHaveCount(1);
});

test('Clear reference controls are available in Single, Joint and Multilayer',async({page})=>{
 await page.goto(base);
 await expect(page.locator('#reference_ri_clear')).toBeVisible();
 await page.locator('#workflow').selectOption('joint');
 await expect(page.locator('#reference_ri_clear')).toBeVisible();
 await page.locator('#workflow').selectOption('multilayer');
 const frame=page.frameLocator('#multilayer_frame');
 await expect(frame.locator('#referenceClear')).toBeAttached();
 await frame.locator('#reference').setInputFiles('examples/v2/transparent_film_120nm.csv');
 await frame.locator('#referenceClear').click();
 await expect(frame.locator('#reference')).toHaveJSProperty('value','');
});

test('Tutorial fit reliability is consolidated and includes calibrated scientific criteria',async({page})=>{
 await page.goto(base);
 await page.locator('#nav_tutorial').click();
 await expect(page.locator('#tutorial-diagnostics')).toBeVisible();
 await expect(page.locator('#tutorial-diagnostics')).toContainText('Jacobian rank');
 await expect(page.locator('#tutorial-diagnostics')).toContainText('Multistart');
 await expect(page.locator('#tutorial-diagnostics')).toContainText('condition-number');
 await expect(page.getByText('Decision guide · Is my fit good enough?')).toHaveCount(0);
});
test('Assisted fit assessment exists across Single, Joint, and embedded Multilayer',async({page})=>{
 await page.goto(base);
 await page.locator('#v2-science-validation').locator('summary').click();
 await page.locator('#v2-validation-tabs button[data-tab="assisted"]').click();
 await expect(page.locator('#v2-validation-content')).toContainText('NOT ASSESSED');
 await page.locator('#workflow').selectOption('joint');
 await expect(page.locator('#v2-validation-tabs button[data-tab="assisted"]')).toHaveCount(1);
 await page.locator('#workflow').selectOption('multilayer');
 const frame=page.frameLocator('#multilayer_frame');
 await frame.locator('#v2-science-validation summary').click();
 await frame.locator('#v2-validation-tabs button[data-tab="assisted"]').click();
 await expect(frame.locator('#v2-validation-content')).toContainText('NOT ASSESSED');
});

test('Assisted assessment never certifies an underidentified but low-RMSE fit',async({page})=>{
 await page.goto(base);
 await page.locator('#v2-science-validation summary').click();
 await page.locator('#v2-validation-tabs button[data-tab="assisted"]').click();
 await page.evaluate(()=>{
  window.oclValidationResult=()=>({
   success:true,weighted_rmse:1e-10,
   sensitivity:{rank:3,parameters:5,condition_number:200},
   wavelength_nm:[500,600],n:[1.5,1.5],k:[0,0],
   predictions:[{R_fit:[0.1,0.12],T_fit:[0.9,0.88]}]
  });
  document.dispatchEvent(new Event('ocl-fit-complete'));
 });
 const report=page.locator('#v2-validation-content');
 await expect(report).toContainText('FAIL (local identifiability)');
 await expect(report).toContainText('3/5 free parameter directions');
 await expect(report).toContainText('Material accuracy');
 await expect(report).toContainText('NOT ASSESSED');
});

test('Joint running status identifies selected Joint model, not stale Single model',async({page})=>{
 await page.goto(base);
 await page.locator('#workflow').selectOption('joint');
 await page.locator('#files').setInputFiles(['examples/v2/transparent_film_120nm.csv','examples/v2/transparent_film_240nm.csv']);
 await page.locator('#joint_method').selectOption('tauc_lorentz');
 await page.evaluate(()=>{window.confirm=()=>true;});
 await page.locator('#run').click();
 await expect(page.locator('#status')).toContainText('Running tauc_lorentz',{timeout:8000});
});
