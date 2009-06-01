import sys
import os
import popen2
import re
import select
import string

comment = re.compile(r'\s*#.*')
assignment = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)')
continuation = re.compile(r'^\s')

INTERSECT = sys.argv[1]
TESTS = sys.argv[2]
TGT_DIR = sys.argv[3]

if not os.path.exists(TGT_DIR):
  os.makedirs(TGT_DIR)

pre = ''
assignments = {}

SHARK = os.path.exists('/usr/bin/shark') and '/usr/bin/shark' or None

def run(cmd):
  test_name, args, op = cmd.split('|')
  test_name = test_name.strip()
  args = args.strip().split()
  op = string.Template(op.strip()).substitute(**assignments)
  print >>sys.stderr, test_name, '...',

  if SHARK:
    CMD = (SHARK, '-o', os.path.join(TGT_DIR, 'prof_%s' % (test_name,)), '-G', '-i', '-1', '-c', '13') + (INTERSECT,) + tuple(args) + (op,)
  else:
    CMD = ('/usr/bin/time', '-v', '-o', os.path.join(TGT_DIR, 'time_%s' % (test_name,))) + (INTERSECT,) + tuple(args) + (op,)

  cmd = popen2.Popen3(CMD, capturestderr = True)
  cout, cerr = cmd.fromchild, cmd.childerr
  cmd.tochild.close()

  r_out = []
  r_err = []
  streams = [ cout, cerr ]
  while len(streams):
    ready = select.select(streams, [], [])

    if cout in ready[0]:
      r_out.append(os.read(cout.fileno(), 1024))
      if not len(r_out[-1]): streams.remove(cout)

    if cerr in ready[0]:
      r_err.append(os.read(cerr.fileno(), 1024))
      if not len(r_err[-1]): streams.remove(cerr)

  exitcode = cmd.wait()

  r_out = ''.join(r_out)
  r_err = ''.join(r_err)

  open (os.path.join(TGT_DIR, 'test_%s.out' % (test_name,)), 'w').write(r_out)
  open (os.path.join(TGT_DIR, 'test_%s.err' % (test_name,)), 'w').write(r_err)
  if exitcode:
    print >>sys.stderr, 'FAIL'
  else:
    print >>sys.stderr, 'PASS'

c = []
def consume():
  if len(c):
    cmd = '\n'.join(c)
    m = assignment.match(cmd)
    if m is not None:
      k, v = m.groups()
      assignments[k] = eval(v)
    else:
      run(re.sub(r'\s+', ' ', cmd))
    c[:] = []

for cmd in [ _ for _ in open(TESTS) if comment.sub('', _).strip() ]:
  cmd = cmd.rstrip('\n')

  if cmd == '':
    consume()
    continue

  if continuation.match(cmd):
    c.append(cmd.lstrip())
    continue

  consume()

  c.append(cmd)

consume()
