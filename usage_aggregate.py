"""Privacy-minimal daily aggregate statistics.

Only (UTC day, event type, count) is persisted. No IP, session ID, cookie,
browser fingerprint, filename, spectral data, or model parameter enters the DB.
"""
from datetime import datetime, timedelta, timezone
import os
from flask import Blueprint, jsonify, request

usage_api = Blueprint('ocl_usage', __name__)
EVENTS = ('visit','fit_single','fit_joint','fit_multilayer','export_csv','export_excel','project_save')
PUBLIC_EVENTS = frozenset(('visit','export_csv','project_save'))
FIT_ENDPOINTS = {
    '/api/original-inverse':'fit_single',
    '/api/classical-analysis':'fit_single',
    '/api/custom-fit':'fit_single',
    '/api/fit':'fit_single',
    '/api/joint-fit':'fit_joint',
    '/api/multilayer/joint-fit':'fit_multilayer',
    '/api/export-causal-workbook':'export_excel',
}
CREATE_SQL = """
CREATE TABLE IF NOT EXISTS ocl_usage_daily (
  day DATE NOT NULL,
  event VARCHAR(24) NOT NULL,
  count BIGINT NOT NULL DEFAULT 0 CHECK (count >= 0),
  PRIMARY KEY (day, event)
)
"""

def configured():
    return bool(os.environ.get('USAGE_DATABASE_URL','').strip())

def connect():
    import psycopg
    return psycopg.connect(os.environ['USAGE_DATABASE_URL'],connect_timeout=4)

def record(event):
    if not configured() or event not in EVENTS:return False
    try:
        with connect() as conn:
            with conn.cursor() as cursor:
                cursor.execute(CREATE_SQL)
                cursor.execute(
                    'INSERT INTO ocl_usage_daily(day,event,count) VALUES (%s,%s,1) '
                    'ON CONFLICT (day,event) DO UPDATE SET count=ocl_usage_daily.count+1',
                    (datetime.now(timezone.utc).date(),event))
        return True
    except Exception:
        # Analytics must never affect inversion or file exports.
        return False

@usage_api.after_app_request
def count_completed_requests(response):
    event=FIT_ENDPOINTS.get(request.path)
    if request.method=='POST' and event and 200<=response.status_code<300:
        record(event)
    return response

@usage_api.post('/api/usage/event')
def usage_event():
    if not request.is_json:return jsonify({'error':'JSON required'}),415
    payload=request.get_json(silent=True)
    if not isinstance(payload,dict) or set(payload)!={'event'} or payload.get('event') not in PUBLIC_EVENTS:
        return jsonify({'error':'Unknown aggregate event'}),400
    if not configured():return jsonify({'recorded':False,'reason':'Persistent storage not configured'}),503
    if not record(payload['event']):return jsonify({'recorded':False,'reason':'Database unavailable'}),503
    return jsonify({'recorded':True})

@usage_api.get('/api/usage/summary')
def usage_summary():
    try:days=int(request.args.get('days','30'))
    except ValueError:return jsonify({'error':'Invalid days'}),400
    if days not in (7,30,365):return jsonify({'error':'Allowed periods: 7, 30, 365 days'}),400
    if not configured():
        return jsonify({'enabled':False,'reason':'Persistent aggregate storage is not connected. No usage history has been collected.','days':days})
    end=datetime.now(timezone.utc).date()
    start=end-timedelta(days=days-1)
    try:
        with connect() as conn:
            with conn.cursor() as cursor:
                cursor.execute(CREATE_SQL)
                cursor.execute('SELECT day,event,count FROM ocl_usage_daily WHERE day >= %s AND day <= %s',(start,end))
                rows=cursor.fetchall()
    except Exception:
        return jsonify({'enabled':False,'reason':'Aggregate database temporarily unavailable','days':days})
    byday={str(start+timedelta(days=i)):{'date':str(start+timedelta(days=i)),**{event:0 for event in EVENTS}} for i in range(days)}
    totals={event:0 for event in EVENTS}
    for day,event,count in rows:
        key=str(day)
        if key in byday and event in totals:
            byday[key][event]+=int(count)
            totals[event]+=int(count)
    response=jsonify({'enabled':True,'days':days,'daily':list(byday.values()),
        'totals':totals,'measurement':'page visits and aggregate actions, not unique people',
        'privacy':'Only day, event category and count are stored'})
    response.headers['Cache-Control']='no-store'
    return response
