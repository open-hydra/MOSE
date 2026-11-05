#!/usr/bin/env python3
"""Per-record and per-cell copies of the multirate integrator (src/lib/numerics/time/explicit/Mod_Multirate.f90) against
the routines of the single-rate stage they copy, statement by statement (with one level the stage
calls the block routines themselves, so the one-level identity test does not cover the copies; this test does):
  SA_Source_Cell   = body of the cell loop of SA_Source_Blk (Lib_Spalart.f90) and the constants of that module;
  SST_Source_Cell  = body of the cell loop of SST_Blk (Lib_SST.f90) and the constants of that module;
  Ghost_Record      = the select of the loop on the local records of Fill_Ghost_Cell (Lib_Ghost.f90): the same case labels
                      calling the same ghost routines (Ghost_Wall at the walls and the 400 faces: one set of ghost values
                      for the fluxes and the sources), a call under a one-line if with its condition, the chimera aside
                      (rejected by the multirate);
  Connection_Record = the body of the loop of Connection_Fluxes (Mod_BC_Fluxes.f90), one flux per connection face;
  BC_Flux_Record    = the branches of the select of BC_Fluxes (Mod_BC_Fluxes.f90) for the boundary types that Check_Domain
                      accepts, statement by statement: domain % bc(i) written bc, the wall branches (301, 302) in their
                      non-rotating form (the rotating frame is rejected), a block if of one statement written on one line;
                      and every type that Check_Domain accepts has its branch (0: no flux, as in BC_Fluxes);
  Dtcap_Cell        = the statements of the internal procedure compute of Compute_dt (Mod_dt.f90) without the branch of
                      the preconditioning (rejected by the multirate), and the wedge lines (rc, er) and the bound of the
                      stiffness of the explicit turbulence sources (dstiff, its flag explicit_sources written out) of
                      the cell loop of Compute_dt, with the arrays of the cell for the arguments
                      (rans_ = blk%P(nt:nprim,i,j,k));
  Check_And_Fix_State = Check_And_Fix_State (Lib_Newstate.f90), the message of the stop marked (multirate);
  Diffusive_Segment = the face metric (the stored one, read in place, or Viscous_Face_Metric), the metrics per group of a
                      face with a wall ghost cell in its stencil (Mg) and the call of Diffusive_Flux of Fluxes_blk
                      (Mod_Fluxes.f90), direction by direction, for the faces in a level and (other residual arrays) for
                      the faces between levels.
And the MPI tags of the multirate (TAG_MR_*) against those of the exchanges of Mod_GhostExchange.f90 (TAG_*): no value
shared.
Statements are compared after removing comments and blank lines and collapsing blanks; continuation lines are joined.
Usage: verify_copies.py SRC_DIR   (exit 0 if every copy is identical to its original)"""
import os
import re
import sys


def statements(lines):
    """Fortran statements of a block of lines: comments and blank lines removed, continuations joined, blanks collapsed,
    lower case outside character constants is not needed (the copies keep the case)."""
    out, cur = [], ''
    for ln in lines:
        s = ln.split('!', 1)[0] if '!' in ln and not ln.lstrip().startswith('!$omp') else ln
        s = s.strip()
        if not s or s.startswith('!'):
            continue
        cont = s.endswith('&')
        if cont:
            s = s[:-1].rstrip()
        if s.startswith('&'):
            s = s[1:].lstrip()
        cur = (cur + ' ' + s).strip()
        if not cont:
            out.append(re.sub(r'\s+', ' ', cur))
            cur = ''
    return out


def block(lines, start_pat, end_pat, after=0):
    a = next(n for n in range(after, len(lines)) if re.search(start_pat, lines[n]))
    b = next(n for n in range(a + 1, len(lines)) if re.search(end_pat, lines[n]))
    return a, b


def loop_body(path, routine):
    L = open(path, encoding='utf-8', errors='replace').read().split('\n')
    h, _ = block(L, r'^\s*subroutine\s+' + routine + r'\b', r'^\s*end subroutine\s+' + routine)
    a, b = block(L, r'^\s*do i = 1, n\(1\)\s*$', r'^\s*enddo ; enddo ; enddo', h)
    return statements(L[a + 1:b])


def constants(path, last):
    L = open(path, encoding='utf-8', errors='replace').read().split('\n')
    a, b = block(L, r'real\(R8\), parameter, private :: &', last)
    return statements(L[a:b + 1])


def copy_body(path, routine, first_executable):
    """statements of the copy after its declarations: from the first executable statement to the end of the routine"""
    L = open(path, encoding='utf-8', errors='replace').read().split('\n')
    h, e = block(L, r'^\s*subroutine\s+' + routine + r'\b', r'^\s*end subroutine\s+' + routine)
    st = statements(L[h + 1:e])
    k = next(n for n, s in enumerate(st) if s.startswith(first_executable))
    return st[:k], st[k:]


def main():
    src = sys.argv[1]
    mr = os.path.join(src, 'lib', 'numerics', 'time', 'explicit', 'Mod_Multirate.f90')
    ok = True
    for routine, orig, blk, last, first in (
            ('SA_Source_Cell', 'lib/physics/turbulence/Lib_Spalart.f90', 'SA_Source_Blk', r'tolerance = 1.d-2',
             'call Co_rotot_Rtot'),
            ('SST_Source_Cell', 'lib/physics/turbulence/Lib_SST.f90', 'SST_Blk', r'cr3 = 1.d0', 'call Co_rotot_Rtot')):
        o = loop_body(os.path.join(src, orig), blk)
        decl, c = copy_body(mr, routine, first)
        same = o == c
        print('  %s: %d statements, body of the cell loop of %s: %s' % (routine, len(c), blk,
                                                                        'identical' if same else 'DIFFERENT'))
        if not same:
            for n, (x, y) in enumerate(zip(o, c)):
                if x != y:
                    print('    first difference at statement %d:\n      original: %s\n      copy:     %s' % (n + 1, x, y))
                    break
            if len(o) != len(c):
                print('    statements: original %d, copy %d' % (len(o), len(c)))
        k = constants(os.path.join(src, orig), last)
        kc = [s for s in decl if s.startswith('real(R8), parameter ::')]
        kc_full = kc[0] if kc else ''
        ko = k[0].replace('real(R8), parameter, private ::', 'real(R8), parameter ::') if k else '?'
        csame = ko == kc_full
        print('  %s: constants of %s: %s' % (routine, os.path.basename(orig), 'identical' if csame else 'DIFFERENT'))
        ok = ok and same and csame
    ok = Ghost_Select(src, mr) and ok
    # Connection_Record against the body of the loop of Connection_Fluxes (one flux per connection face), after its
    # three selections (record index, connection type, owner record)
    B = open(os.path.join(src, 'lib', 'numerics', 'fluxes', 'bc', 'Mod_BC_Fluxes.f90'), encoding='utf-8',
             errors='replace').read().split('\n')
    h, _ = block(B, r'^\s*subroutine\s+Connection_Fluxes\b', r'^\s*end subroutine\s+Connection_Fluxes')
    a, b = block(B, r'^\s*do ii = 1, domain % n_local_bc\s*$', r'^\s*enddo\s*$', h)
    co = statements(B[a + 1:b])[3:]
    _, cc = copy_body(mr, 'Connection_Record', 'Bm = domain % bc(i) % b')
    csame = co == cc
    print('  Connection_Record: %d statements, body of the loop of Connection_Fluxes: %s'
          % (len(cc), 'identical' if csame else 'DIFFERENT'))
    if not csame:
        for n, (x, y) in enumerate(zip(co, cc)):
            if x != y:
                print('    first difference at statement %d:\n      original: %s\n      copy:     %s' % (n + 1, x, y))
                break
    ok = ok and csame
    ok = BC_Select(src, mr) and ok
    ok = Dtcap_Copy(src, mr) and ok
    ok = Check_State_Copy(src, mr) and ok
    ok = Diffusive_Calls(src, mr) and ok
    ok = Tags(src, mr) and ok
    print('verify_copies: %s' % ('ok' if ok else 'FAIL'))
    return 0 if ok else 1


def first_difference(o, c):
    for n, (x, y) in enumerate(zip(o, c)):
        if x != y:
            print('    first difference at statement %d:\n      original: %s\n      copy:     %s' % (n + 1, x, y))
            return
    if len(o) != len(c):
        print('    statements: original %d, copy %d' % (len(o), len(c)))


def select_cases(st, start):
    """case label (blanks removed) -> statements of its branch, from the select whose statement starts with start to its
    end select (no select nested in it)"""
    k = next(n for n, x in enumerate(st) if x.startswith(start))
    d, lab = {}, None
    for x in st[k + 1:]:
        if x.startswith('end select'):
            break
        if x.startswith('case'):
            lab = x.replace(' ', '')
            d[lab] = []
        elif lab is not None:
            d[lab].append(x)
    return d


def case_types(lab):
    """the boundary types of a case label such as case(0,101,201,300,301,302,400:408,410,420)"""
    out = set()
    for item in lab[len('case('):-1].split(','):
        if ':' in item:
            a, b = item.split(':')
            out |= set(range(int(a), int(b) + 1))
        else:
            out.add(int(item))
    return out


def BC_Select(src, mr):
    """BC_Flux_Record against the select of BC_Fluxes (Mod_BC_Fluxes.f90): for every case of the copy, the branch of the
    same case of BC_Fluxes with domain % bc(i) written bc, the branch of the rotating frame of the walls removed (if ... then
    the rotating form, elseif (model>0) then the other: the other, under if (model>0)), a block if of one statement
    written on one line; blanks ignored. The types that Check_Domain accepts (the cases without an unconditional Reject)
    must each have a case in the copy."""
    def norm(st):
        out, k = [], 0
        st = [x.replace('domain % bc(i)', 'bc') for x in st]
        while k < len(st):
            x = st[k]
            if x.startswith('if (model>0.AND.obj_rot%enabled'):
                e = next(n for n in range(k + 1, len(st)) if st[n] == 'elseif (model>0) then')
                f = next(n for n in range(e + 1, len(st)) if st[n] == 'endif')
                out += ['if (model>0) ' + y for y in st[e + 1:f]]
                k = f + 1
                continue
            if x.startswith('if ') and x.endswith(' then') and k + 2 < len(st) and st[k + 2] in ('endif', 'end if'):
                out.append(x[:-len(' then')] + ' ' + st[k + 1])
                k += 3
                continue
            out.append(x)
            k += 1
        return [x.replace(' ', '') for x in out]
    B = open(os.path.join(src, 'lib', 'numerics', 'fluxes', 'bc', 'Mod_BC_Fluxes.f90'), encoding='utf-8',
             errors='replace').read().split('\n')
    h, e = block(B, r'^\s*subroutine\s+BC_Fluxes\b', r'^\s*end subroutine\s+BC_Fluxes')
    o = select_cases(statements(B[h + 1:e]), 'select case ( domain % bc(i) % type )')
    L = open(mr, encoding='utf-8', errors='replace').read().split('\n')
    h, e = block(L, r'^\s*subroutine\s+BC_Flux_Record\b', r'^\s*end subroutine\s+BC_Flux_Record')
    c = select_cases(statements(L[h + 1:e]), 'select case ( bc % type )')
    h, e = block(L, r'^\s*subroutine\s+Check_Domain\b', r'^\s*end subroutine\s+Check_Domain')
    cd = select_cases(statements(L[h + 1:e]), 'select case (t)')
    accepted = set()
    for lab, st in cd.items():
        if lab != 'casedefault' and not (st and st[0].startswith('call Reject')):
            accepted |= case_types(lab)
    ok, n = True, 0
    covered = set()
    for lab, st in c.items():
        if lab == 'casedefault':
            continue
        covered |= case_types(lab)
        if lab == 'case(0)':
            if st:
                print('    %s: the copy has statements, BC_Fluxes gives these types no flux' % lab)
                ok = False
            continue
        n += 1
        if lab not in o:
            print('    %s: no such case in BC_Fluxes' % lab)
            ok = False
            continue
        a, b = norm(o[lab]), norm(st)
        if a != b:
            print('    %s: DIFFERENT' % lab)
            first_difference(a, b)
            ok = False
    missing = sorted(accepted - covered)
    if missing:
        print('    types accepted by Check_Domain without a case in BC_Flux_Record: %s' % missing)
        ok = False
    print('  BC_Flux_Record: %d cases, the select of BC_Fluxes for the types accepted by Check_Domain (%d): %s'
          % (n, len(accepted), 'identical' if ok else 'DIFFERENT'))
    return ok


def Dtcap_Copy(src, mr):
    """Dtcap_Cell against compute of Compute_dt and the wedge lines of its cell loop (Mod_dt.f90)."""
    D = open(os.path.join(src, 'lib', 'numerics', 'time', 'Mod_dt.f90'), encoding='utf-8', errors='replace').read().split('\n')
    h, e = block(D, r'^\s*subroutine\s+compute\b', r'^\s*end subroutine\s+compute')
    st = statements(D[h + 1:e])
    st = st[next(n for n, x in enumerate(st) if x.startswith('call co_rotot_Rtot')):]
    o, skip = [], 0
    for x in st:                                     # the branch of the preconditioning out, its else branch in
        if x == 'if ( obj_prec%enabled ) then':
            skip = 1
            continue
        if skip == 1 and x == 'else':
            skip = 2
            continue
        if skip and x == 'endif':
            skip = 0
            continue
        if skip == 1:
            continue
        o.append(x)
    hc, ec = block(D, r'^\s*subroutine\s+Compute_dt\b', r'^\s*end subroutine\s+Compute_dt')
    cl = statements(D[hc + 1:ec])
    w = [x for x in cl if x.startswith('wedge = ')]
    a = next(n for n, x in enumerate(cl) if x == 'rc = 0d0')
    b = next(n for n, x in enumerate(cl) if x == 'end if' and n > a)
    o = w + [x.replace('domain%blk(b)%node(', 'blk%node(') for x in cl[a:b + 1]] + o
    o = [x.replace('rans_variables=rans_', 'rans_variables=blk%P(nt:nprim,i,j,k)').replace('rans_(1)', 'blk%P(nt,i,j,k)')
         for x in o]
    L = open(mr, encoding='utf-8', errors='replace').read().split('\n')
    h, e = block(L, r'^\s*real\(R8\) function\s+Dtcap_Cell\b', r'^\s*end function\s+Dtcap_Cell')
    c = statements(L[h + 1:e])
    # the bound of the stiffness of the explicit turbulence sources: the two statements of the cell loop of Compute_dt
    # (dstiff), with its flag explicit_sources written out and the arrays of the cell
    flag = [x for x in cl if x.startswith('explicit_sources = ')]
    a = next((n for n, x in enumerate(cl) if x == 'dstiff = 0d0'), None)
    so = []
    if len(flag) == 1 and a is not None and a + 1 < len(cl):
        cond = flag[0].split('=', 1)[1].strip()
        so = [x.replace('if ( explicit_sources )', 'if ( %s )' % cond).replace('dstiff', 'stiff')
              .replace('domain%blk(b)%', 'blk%') for x in cl[a:a + 2]]
    a = next((n for n, x in enumerate(c) if x == 'stiff = 0d0'), None)
    sc = c[a:a + 2] if a is not None else []
    ssame = len(so) == 2 and so == sc
    print('  Dtcap_Cell: the bound of the stiffness of the turbulence sources of the cell loop of Compute_dt: %s'
          % ('identical' if ssame else 'DIFFERENT'))
    if not ssame:
        first_difference(so, sc)
    k = [n for n, x in enumerate(c) if x.startswith('wedge = ') or x.startswith('call co_rotot_Rtot')]
    c = c[k[0]:] if k else c
    same = o == c
    print('  Dtcap_Cell: %d statements, compute and the wedge lines of Compute_dt: %s' % (len(c), 'identical' if same
                                                                                           else 'DIFFERENT'))
    if not same:
        first_difference(o, c)
    return same and ssame


def Check_State_Copy(src, mr):
    """Check_And_Fix_State of the multirate against the one of Lib_Newstate.f90 (the message of the stop aside)."""
    N = open(os.path.join(src, 'lib', 'numerics', 'time', 'explicit', 'Lib_Newstate.f90'), encoding='utf-8',
             errors='replace').read().split('\n')
    h, e = block(N, r'^\s*subroutine\s+Check_And_Fix_State\b', r'^\s*end subroutine\s+Check_And_Fix_State')
    st = statements(N[h + 1:e])
    o = st[next(n for n, x in enumerate(st) if x.startswith('if ( any(')):]
    _, c = copy_body(mr, 'Check_And_Fix_State', 'if ( any(')
    c = [x.replace(' (multirate)"', '"') for x in c]
    same = o == c
    print('  Check_And_Fix_State: %d statements, Check_And_Fix_State of Lib_Newstate: %s' % (len(c), 'identical' if same
                                                                                           else 'DIFFERENT'))
    if not same:
        first_difference(o, c)
    return same


def Diffusive_Calls(src, mr):
    """Face metric and call of Diffusive_Flux in Diffusive_Segment against Fluxes_blk, per direction: the call for a face
    in the level is the one of Fluxes_blk; the call for a face between levels differs only in the residual arrays."""
    F = open(os.path.join(src, 'lib', 'numerics', 'fluxes', 'Mod_Fluxes.f90'), encoding='utf-8', errors='replace').read()
    F = F.split('\n')
    h, e = block(F, r'^\s*subroutine\s+Fluxes_blk\b', r'^\s*end subroutine\s+Fluxes_blk')
    o = Face_Blocks(statements(F[h + 1:e]))
    _, c = copy_body(mr, 'Diffusive_Segment', 'Sc = ')
    c = Face_Blocks(c)
    same = len(o) == 3 and len(c) == 6
    for d in range(3) if same else ():
        ob = o[d]
        cb1, cb2 = c[2 * d:2 * d + 2]
        res = re.search(r'blk % R\(:,[^)]*\), blk % R\(:,[^)]*\)', ob[-1])
        same = same and res is not None and ob == cb1 and \
            ob[:-1] + [ob[-1].replace(res.group(0), 'Rt(:,0), Rt(:,1)')] == cb2
    print('  Diffusive_Segment: %d calls, face metrics and call of Diffusive_Flux of Fluxes_blk: %s'
          % (len(c), 'identical' if same else 'DIFFERENT'))
    return same


def Face_Blocks(st):
    """the statements of each viscous face: from the face metric (nullify ( Mg ), then the stored metric of the face or
    Mf = Viscous_Face_Metric) to the call of Diffusive_Flux"""
    out, cur = [], None
    for x in st:
        if cur is None and (x.startswith('nullify ( Mg )') or x.startswith('Mf = Viscous_Face_Metric')):
            cur = [x]
        elif cur is not None:
            cur.append(x)
            if x.startswith('call Diffusive_Flux'):
                out.append(cur)
                cur = None
    return out


def Ghost_Select(src, mr):
    """the case labels of the select of Ghost_Record and the ghost routines each case calls, against the select of the
    loop on the local records of Fill_Ghost_Cell (Lib_Ghost.f90); the chimera (102) aside"""
    def cases(st):
        """case label -> the routines its branch calls, a call under a one-line if with its condition"""
        d, lab = {}, None
        for x in st:
            if x.startswith('select case'):
                continue
            if x.startswith('end select'):
                break
            if x.startswith('case'):
                lab = x.replace(' ', '')
                d[lab] = []
            elif lab is not None and x.startswith('call '):
                d[lab].append(x.split()[1])
            elif lab is not None and x.startswith('if') and ' call ' in x:
                cond, rest = x.split(' call ', 1)                   # a call under a one-line if: the condition too
                d[lab].append(cond.replace(' ', '') + ' ' + rest.split()[0])
        d.pop('case(102)', None)
        return d
    G = open(os.path.join(src, 'lib', 'numerics', 'space', 'Lib_Ghost.f90'), encoding='utf-8', errors='replace').read()
    G = G.split('\n')
    h, e = block(G, r'^\s*subroutine\s+Fill_Ghost_Cell\b', r'^\s*end subroutine\s+Fill_Ghost_Cell')
    st = statements(G[h + 1:e])
    k = next(n for n, x in enumerate(st) if x.startswith('select case ( domain % bc(i) % type'))
    o = cases(st[k:])
    _, c = copy_body(mr, 'Ghost_Record', 'Bm = ')
    k = next(n for n, x in enumerate(c) if x.startswith('select case'))
    c = cases(c[k:])
    same = o == c and len(o) > 0
    print('  Ghost_Record: %d cases, the select of the local records of Fill_Ghost_Cell: %s'
          % (len(c), 'identical' if same else 'DIFFERENT'))
    if not same:
        for lab in sorted(set(o) | set(c)):
            if o.get(lab) != c.get(lab):
                print('    %s: original %s, copy %s' % (lab, o.get(lab), c.get(lab)))
    return same


def Tags(src, mr):
    """the values of the MPI tags of the multirate (integer parameters TAG_MR_* of Mod_Multirate.f90) against those of the
    exchanges of the single-rate code (TAG_* of Mod_GhostExchange.f90): no value in common, so that a receive of one
    cannot take a message of the other"""
    def tags(path, prefix):
        t = {}
        for s in statements(open(path, encoding='utf-8', errors='replace').read().split('\n')):
            if s.startswith('integer, parameter') and '::' in s:
                for item in s.split('::', 1)[1].split(','):
                    if '=' in item:
                        name, val = (x.strip() for x in item.split('=', 1))
                        if name.startswith(prefix) and re.fullmatch(r'-?\d+', val):
                            t[name] = int(val)
        return t
    g = tags(os.path.join(src, 'lib', 'parallel', 'Mod_GhostExchange.f90'), 'TAG_')
    m = tags(mr, 'TAG_MR_')
    shared = sorted(set(g.values()) & set(m.values()))
    ok = len(g) > 0 and len(m) > 0 and not shared
    print('  MPI tags: multirate %s, Mod_GhostExchange %s: %s' % (
        ' '.join(str(v) for v in sorted(m.values())), ' '.join(str(v) for v in sorted(g.values())),
        'no value shared' if ok else ('SHARED ' + ' '.join(str(v) for v in shared) if shared else 'NOT FOUND')))
    return ok


if __name__ == '__main__':
    sys.exit(main())
