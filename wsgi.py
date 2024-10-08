
python_home='/var/www/DESIRE/venv/desire_py36'
activate_this=python_home+'/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))

def application(environ, start_response):
    status = '200 OK'
    output = b'Hello World!'

    response_headers = [('Content-type', 'text/plain'),
                        ('Content-Length', str(len(output)))]
    start_response(status, response_headers)

    return [output]
