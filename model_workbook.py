"""Formula-based Excel reconstruction of fitted causal dielectric models."""
from io import BytesIO
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter

def build_model_workbook(model, coefficients, wavelengths_nm, metadata=None):
    if model not in ('lorentz','drude_lorentz'):
        raise ValueError('Excel model reconstruction requires Lorentz or Drude-Lorentz coefficients')
    offset=3 if model=='drude_lorentz' else 1
    count=(len(coefficients)-offset)//3
    if not (1<=count<=8) or offset+3*count!=len(coefficients):
        raise ValueError('Invalid causal parameter count')
    wb=Workbook()
    p=wb.active;p.title='Parameters'
    p.append(['Parameter','Value','Units'])
    p.append(['epsilon_inf',float(coefficients[0]),'dimensionless'])
    if model=='drude_lorentz':
        p.append(['Drude plasma energy',float(coefficients[1]),'eV'])
        p.append(['Drude damping',float(coefficients[2]),'eV'])
    for j in range(count):
        for name,value,unit in zip((f'osc{j+1}_strength',f'osc{j+1}_energy',f'osc{j+1}_damping'),coefficients[offset+3*j:offset+3*j+3],('dimensionless','eV','eV')):
            p.append([name,float(value),unit])
    sheet=wb.create_sheet('Reconstruction')
    sheet.append(['wavelength_nm','photon_energy_eV','epsilon_real','epsilon_imag','n','k'])
    for row,w in enumerate(wavelengths_nm,2):
        sheet.cell(row,1,float(w))
        sheet.cell(row,2,f'=1239.841984/A{row}')
        real_terms=[];imag_terms=[]
        for j in range(count):
            first=offset+3*j+2 # Excel parameter row (zero-based coefficient index + 2)
            f=f'Parameters!$B${first}'
            e=f'Parameters!$B${first+1}'.replace('\\$','$')
            g=f'Parameters!$B${first+2}'.replace('\\$','$')
            a=f'({e}^2-B{row}^2)'
            b=f'({g}*B{row})'
            denom=f'({a}^2+{b}^2)'
            real_terms.append(f'({f}*{e}^2*{a}/{denom})')
            imag_terms.append(f'({f}*{e}^2*{b}/{denom})')
        if model=='drude_lorentz':
            wp='Parameters!$B$3';gd='Parameters!$B$4'
            real_terms.append(f'(-{wp}^2/(B{row}^2+{gd}^2))')
            imag_terms.append(f'({wp}^2*{gd}/(B{row}*(B{row}^2+{gd}^2)))')
        sheet.cell(row,3,'=Parameters!$B$2'+''.join('+'+t for t in real_terms))
        sheet.cell(row,4,'='+('+'.join(imag_terms) if imag_terms else '0'))
        sheet.cell(row,5,f'=SQRT((SQRT(C{row}^2+D{row}^2)+C{row})/2)')
        sheet.cell(row,6,f'=SQRT(MAX(0,(SQRT(C{row}^2+D{row}^2)-C{row})/2))')
    meta=wb.create_sheet('Metadata')
    meta.append(['Field','Value'])
    meta.append(['Model',model]);meta.append(['Oscillators',count])
    for k,v in (metadata or {}).items():meta.append([str(k),str(v)])
    notes=wb.create_sheet('Instructions')
    for line in [
        'Edit wavelength_nm values in Reconstruction; Excel recalculates epsilon, n and k.',
        'E[eV] = 1239.841984 / wavelength[nm].',
        'Lorentz epsilon = epsilon_inf + SUM(f*E0^2/(E0^2-E^2-i*gamma*E)).',
        'Drude addition = -plasma_energy^2/(E^2+i*gamma_D*E).',
        'Convention: complex refractive index n+i*k, epsilon=(n+i*k)^2.',
        'Formula results may not be cached until opened in Excel or LibreOffice.',
        'Outside the fitted spectral window, reconstruction is extrapolation, not validated measurement.'
    ]:notes.append([line])
    for ws in wb:
        ws.freeze_panes='A2'
        for cell in ws[1]:
            cell.font=Font(bold=True,color='FFFFFF')
            cell.fill=PatternFill('solid',fgColor='17365D')
        for col in ws.columns:
            letter=get_column_letter(col[0].column)
            ws.column_dimensions[letter].width=min(88,max(16,max(len(str(c.value or '')) if c.data_type!='f' else 20 for c in col)+3))
    output=BytesIO();wb.save(output);output.seek(0);return output


def build_numerical_model_workbook(model, coefficients, wavelengths_nm, metadata=None):
    """Portable, explicitly numerical model reconstruction for nonlinear dispersion families.

    Coefficients and executable Python instructions are supplied.  ModelValues are
    computed from those coefficients and are NOT falsely marked as Excel formulas.
    """
    import numpy as np
    from organic_dispersion import MODEL_PARAMETER_COUNTS, organic_dielectric
    from optical_inversion_package.inversion import spline_design
    from scientific_diagnostics import dielectric
    w=np.asarray(wavelengths_nm,float)
    p=np.asarray(coefficients,float)
    if w.ndim!=1 or p.ndim!=1 or not len(w) or not np.all(np.isfinite(w)) or np.any(w<=0) or not np.all(np.isfinite(p)):
        raise ValueError('Invalid wavelengths or model coefficients')
    meta=dict(metadata or {})
    if model in MODEL_PARAMETER_COUNTS:
        count=MODEL_PARAMETER_COUNTS[model]
        if len(p)<1 or (len(p)-1)%count:raise ValueError('Incorrect oscillator parameters')
        oscillators=(len(p)-1)//count
        eps=organic_dielectric(w,p,model,oscillators)
        runner='from organic_dispersion import organic_dielectric\\neps=organic_dielectric(wavelength_nm,parameters,model,oscillators)\\nN=np.sqrt(eps+0j)\\nn,k=N.real,N.imag'
    elif model=='kk_constrained':
        centers=np.asarray(meta.get('centers_eV',[]),float)
        widths=np.asarray(meta.get('widths_eV',[]),float)
        if len(centers)!=len(p)-1 or len(widths)!=len(centers):
            raise ValueError('KK basis centers and widths required for parameter reconstruction')
        energy=1239.841984/w
        eps=np.full(w.shape,p[0],complex)
        for strength,center,width in zip(p[1:],centers,widths):
            eps+=strength*center**2/(center**2-energy**2-1j*width*energy)
        runner='energy=1239.841984/wavelength_nm\\neps=np.full(energy.shape,parameters[0],complex)\\nfor strength,center,width in zip(parameters[1:],centers_eV,widths_eV): eps+=strength*center**2/(center**2-energy**2-1j*width*energy)\\nN=np.sqrt(eps+0j)\\nn,k=N.real,N.imag'
    elif model=='spline':
        lo=float(meta['domain_min_nm']);hi=float(meta['domain_max_nm'])
        count=len(p)//2
        if len(p)%2 or count<4:raise ValueError('Spline requires paired n/k coefficients')
        B=spline_design(w,lo,hi,count)
        n=B@p[:count];k=B@p[count:]
        eps=(n+1j*k)**2
        runner='from optical_inversion_package.inversion import spline_design\\nB=spline_design(wavelength_nm,domain_min_nm,domain_max_nm,len(parameters)//2)\\nn=B@parameters[:len(parameters)//2]\\nk=B@parameters[len(parameters)//2:]'
    else:
        raise ValueError('No reproducible numerical reconstruction for this model')
    z=np.sqrt(eps+0j)
    n,k=z.real,z.imag
    wb=Workbook()
    params=wb.active;params.title='Parameters'
    params.append(['Index','Coefficient','Value'])
    for i,v in enumerate(p):params.append([i,'p'+str(i),float(v)])
    vals=wb.create_sheet('Reconstruction')
    vals.append(['wavelength_nm','epsilon_real','epsilon_imag','n','k','alpha_cm_inverse'])
    for wl,e,nv,kv in zip(w,eps,n,k):
        vals.append([float(wl),float(e.real),float(e.imag),float(nv),float(kv),float(4*np.pi*kv/(wl*1e-7))])
    ms=wb.create_sheet('Metadata');ms.append(['Field','Value']);ms.append(['model',model])
    for key,val in meta.items():ms.append([str(key),str(val)])
    doc=wb.create_sheet('External_reconstruction')
    lines=[
        'The Reconstruction sheet is calculated at export; those numbers are not live Excel formula evaluations.',
        'To reconstruct at new wavelengths, run the Python calculation specified below.',
        'Required: Python, numpy, scipy and the matching open-source Optical Constants Laboratory numerical routines.',
        'Wavelength is nm. Parameters in Parameters sheet are in model-native documented units.',
        'This file includes complete numeric coefficients, domain and the algorithm call but does not embed a Python runtime in Excel.',
        'Required call:',
        runner.replace('\\n','; '),
        'Only evaluate inside the fitted wavelength band; extrapolated constants are not validated.',
    ]
    for line in lines:doc.append([line])
    # Include runnable source instructions so changing wavelength does not depend on a Lab server.
    runner_sheet=wb.create_sheet('Python_runner')
    script=[
        'import numpy as np',
        'wavelength_nm=np.array('+repr(w.tolist())+',dtype=float)',
        'parameters=np.array('+repr(p.tolist())+',dtype=float)',
        'model='+repr(model),
        'oscillators='+str((len(p)-1)//MODEL_PARAMETER_COUNTS[model] if model in MODEL_PARAMETER_COUNTS else 0),
    ]
    if model=='kk_constrained':
        script.extend(['centers_eV=np.array('+repr(centers.tolist())+')','widths_eV=np.array('+repr(widths.tolist())+')'])
    if model=='spline':
        script.extend(['domain_min_nm='+repr(lo),'domain_max_nm='+repr(hi),
         'from scipy.interpolate import BSpline',
         'def spline_design(x,low,high,n_basis):',
         '    degree=3',
         '    internal=np.linspace(low,high,n_basis-degree+1)[1:-1]',
         '    knots=np.r_[np.repeat(low,degree+1),internal,np.repeat(high,degree+1)]',
         '    return BSpline(knots,np.eye(n_basis),degree)(x)'])
    script.extend(runner.split('\\n'))
    script.append('for row in zip(wavelength_nm,n,k):print(*row,sep=\",\")')
    for line in script: runner_sheet.append([line])
    if model in MODEL_PARAMETER_COUNTS:
        # The matching material evaluator source is embedded, not linked to a future server.
        from pathlib import Path
        evaluator=wb.create_sheet('organic_dispersion_py')
        for line in Path(__file__).with_name('organic_dispersion.py').read_text(encoding='utf-8').splitlines():
            evaluator.append([line])
        doc.append(['Save organic_dispersion_py column A as organic_dispersion.py beside Python_runner.py.'])
        doc.append(['Install numpy and scipy before executing the two local scripts.'])
    else:
        doc.append(['Save Python_runner column A as Python_runner.py and run with numpy and scipy.'])
    for ws in wb:
        ws.freeze_panes='A2'
        for cell in ws[1]:
            cell.font=Font(bold=True,color='FFFFFF')
            cell.fill=PatternFill('solid',fgColor='17365D')
        for col in ws.columns:
            letter=get_column_letter(col[0].column)
            ws.column_dimensions[letter].width=min(100,max(18,max(len(str(cell.value or '')) for cell in col)+2))
    output=BytesIO();wb.save(output);output.seek(0);return output
