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
