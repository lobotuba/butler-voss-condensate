# -*- coding: utf-8 -*-
"""Build the COMPLETE-PROJECT Butler-Voss Condensate white paper (WP-06) as .docx,
with matplotlib-rendered figures."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

TEAL, RUST, MUT = "#009e8a", "#c2542a", "#61706a"
TEAL_D, RUST_D = RGBColor(0x05, 0x7d, 0x6d), RGBColor(0xac, 0x48, 0x22)
INK, GREY = RGBColor(0x1a, 0x20, 0x1d), RGBColor(0x61, 0x70, 0x6a)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = r"Z:\Butler-Voss Condensate\Butler-Voss Condensate - Whitepaper.docx"
OLD = r"Z:\Butler-Voss Condensate\Butler-Voss Condensate - Screening Whitepaper.docx"

# ----------------------------------------------------------------- figures ----
plt.rcParams.update({"font.family": "serif", "font.serif": ["Georgia", "DejaVu Serif"],
                     "font.size": 11, "axes.edgecolor": "#c5cbc3", "axes.linewidth": .8,
                     "figure.facecolor": "white", "axes.facecolor": "white"})

def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors=MUT, labelsize=9.5, length=3)
    ax.grid(axis="y", color="#e3e6e1", lw=.8); ax.set_axisbelow(True)
    for s in ("left", "bottom"): ax.spines[s].set_color("#c5cbc3")

def fig_drho(p):
    r = [.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5]
    d = [.0537,.0525,.0491,.0441,.0365,.0272,.0196,.0140,.0101,.0077,.0063,.0055]
    fr,fy=[.5,2.5,4.5,6.5,8.5,10.5,11.5],[.0545,.0470,.0345,.0230,.0135,.0075,.0058]
    f,ax=plt.subplots(figsize=(6.4,3.2))
    ax.plot(fr,fy,"--",color=MUT,lw=1.4,label="fit: exp(-r/lambda), lambda = 3.3")
    ax.plot(r,d,"-o",color=RUST,lw=2,ms=5,mfc="white",mec=RUST,mew=1.8,label="measured  drho(r)")
    ax.set_xlabel("radius  r"); ax.set_ylabel("drho(r)"); ax.set_ylim(0,.06); ax.set_xlim(0,12)
    _style(ax); ax.legend(frameon=False,fontsize=9,loc="upper right")
    f.tight_layout(); f.savefig(p,dpi=200,bbox_inches="tight"); plt.close(f)

def fig_box(p):
    x=[10,14,20]
    f,ax=plt.subplots(figsize=(6.4,3.3))
    ax.plot(x,[2.45,3.18,4.24],"-o",color=TEAL,lw=2,ms=6,mfc="white",mec=TEAL,mew=1.8)
    ax.plot(x,[1.86,2.02,2.18],"-o",color=RUST,lw=2,ms=6,mfc="white",mec=RUST,mew=1.8)
    ax.annotate("massless: 1/r (grows)",(20,4.24),xytext=(8,4),textcoords="offset points",
                color=TEAL,fontsize=9.5,fontweight="bold",va="center")
    ax.annotate("massive: 1/m (pinned)",(20,2.18),xytext=(8,0),textcoords="offset points",
                color=RUST,fontsize=9.5,fontweight="bold",va="center")
    ax.set_xlabel("box size"); ax.set_ylabel("apparent range  lambda")
    ax.set_xlim(9,24); ax.set_ylim(1.5,4.6); ax.set_xticks(x); _style(ax)
    f.tight_layout(); f.savefig(p,dpi=200,bbox_inches="tight"); plt.close(f)

def fig_gauge(p):
    e=[.15,.20,.30]; lam=[4.14,3.13,1.94]; gx=[.13,.15,.2,.25,.3,.33]
    f,ax=plt.subplots(figsize=(6.4,3.2))
    ax.plot(gx,[.6/v for v in gx],"--",color=MUT,lw=1.4,label="guide: lambda_L ~ 1/e")
    ax.plot(e,lam,"o",color=RUST,ms=7,mfc="white",mec=RUST,mew=2,label="measured lambda_L")
    ax.set_xlabel("gauge coupling  e"); ax.set_ylabel("penetration depth  lambda_L")
    ax.set_xlim(.12,.33); ax.set_ylim(1.5,4.6); ax.set_xticks(e); _style(ax)
    ax.legend(frameon=False,fontsize=9,loc="upper right")
    f.tight_layout(); f.savefig(p,dpi=200,bbox_inches="tight"); plt.close(f)

def fig_monopole(p):
    r=np.array([2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5])
    B=np.array([.096,.044,.026,.017,.012,.009,.007,.006,.005])
    guide=B[1]*(r[1]/r)**2
    f,ax=plt.subplots(figsize=(6.4,3.2))
    ax.plot(r,guide,"--",color=MUT,lw=1.4,label="reference: 1/r^2")
    ax.plot(r,B,"o",color=TEAL,ms=6,mfc="white",mec=TEAL,mew=1.8,label="monopole field |B|(r)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("radius  r  (log)"); ax.set_ylabel("|B|  (log)")
    _style(ax); ax.grid(which="both",axis="both",color="#e9ece7",lw=.7)
    ax.legend(frameon=False,fontsize=9,loc="lower left")
    f.tight_layout(); f.savefig(p,dpi=200,bbox_inches="tight"); plt.close(f)

def fig_lorentz(p):
    kf=np.array([.05,.1,.2,.3,.4,.6,.8])
    hexa=np.array([1.06e-7,1.70e-6,2.73e-5,1.40e-4,4.52e-4,2.41e-3,8.25e-3])
    fcc=np.array([1.71e-4,6.86e-4,2.75e-3,6.23e-3,1.11e-2,2.56e-2,4.69e-2])
    f,ax=plt.subplots(figsize=(6.4,3.2))
    ax.plot(kf,fcc,"-o",color=RUST,lw=2,ms=5,mfc="white",mec=RUST,mew=1.8,label="fcc 3D  (~ (k/kmax)^2)")
    ax.plot(kf,hexa,"-o",color=TEAL,lw=2,ms=5,mfc="white",mec=TEAL,mew=1.8,label="hex 2D  (~ (k/kmax)^4)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("k / k_max   (energy / E_Planck)"); ax.set_ylabel("anisotropy  dc/c  (Lorentz violation)")
    _style(ax); ax.grid(which="both",axis="both",color="#e9ece7",lw=.7)
    ax.legend(frameon=False,fontsize=9,loc="upper left")
    f.tight_layout(); f.savefig(p,dpi=200,bbox_inches="tight"); plt.close(f)

figs = {}
for nm, fn in [("drho",fig_drho),("box",fig_box),("gauge",fig_gauge),("mono",fig_monopole),
               ("lorentz",fig_lorentz)]:
    figs[nm]=os.path.join(HERE,f"wpf_{nm}.png"); fn(figs[nm])

# ------------------------------------------------------------------- docx -----
doc = Document()
st = doc.styles["Normal"]; st.font.name="Cambria"; st.font.size=Pt(11)
st.paragraph_format.space_after=Pt(7); st.paragraph_format.line_spacing=1.16
for hs, sz in (("Heading 1",15),("Heading 2",12.5),("Title",25)):
    s=doc.styles[hs]; s.font.name="Cambria"; s.font.size=Pt(sz); s.font.color.rgb=INK

def shade(par, hexcol):
    sh=OxmlElement("w:shd"); sh.set(qn("w:val"),"clear"); sh.set(qn("w:fill"),hexcol)
    par._p.get_or_add_pPr().append(sh)

def body(text, justify=True, italic=False, size=None, color=None, after=7, before=0):
    p=doc.add_paragraph(); r=p.add_run(text)
    if italic: r.italic=True
    if size: r.font.size=Pt(size)
    if color: r.font.color.rgb=color
    if justify: p.alignment=WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.space_after=Pt(after); p.paragraph_format.space_before=Pt(before)
    return p

def heading(text, lvl=1):
    h=doc.add_heading(level=lvl); h.paragraph_format.space_before=Pt(12 if lvl==1 else 8)
    parts=text.split("  ",1)
    r=h.add_run(parts[0]+"  "); r.font.color.rgb=TEAL_D
    h.add_run(parts[1] if len(parts)>1 else "")
    return h

def add_eq(markup, num):
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before=Pt(6); p.paragraph_format.space_after=Pt(8)
    i=0
    while i<len(markup):
        c=markup[i]
        if c in "_^" and i+1<len(markup) and markup[i+1]=="{":
            j=markup.index("}",i); r=p.add_run(markup[i+2:j]); r.italic=True; r.font.name="Cambria"
            r.font.subscript = (c=="_"); r.font.superscript = (c=="^"); i=j+1
        else:
            j=i
            while j<len(markup) and not (markup[j] in "_^" and j+1<len(markup) and markup[j+1]=="{"): j+=1
            r=p.add_run(markup[i:j]); r.italic=True; r.font.name="Cambria"; i=j
    t=p.add_run("      ("+num+")"); t.font.color.rgb=GREY
    return p

def caption(no, text):
    p=doc.add_paragraph(); p.paragraph_format.space_after=Pt(12)
    r=p.add_run("Figure "+no+". "); r.bold=True; r.font.size=Pt(9.5)
    r2=p.add_run(text); r2.font.size=Pt(9.5); r2.font.color.rgb=GREY
    p.alignment=WD_ALIGN_PARAGRAPH.JUSTIFY

def add_figure(path, no, cap):
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before=Pt(8); p.paragraph_format.space_after=Pt(4)
    p.add_run().add_picture(path, width=Inches(5.6)); caption(no, cap)

def result(tag, text, warn=False):
    p=doc.add_paragraph(); p.paragraph_format.space_before=Pt(6); p.paragraph_format.space_after=Pt(9)
    p.paragraph_format.left_indent=Inches(0.12); shade(p, "FBEDE6" if warn else "E6F4F1")
    r=p.add_run(tag+"  "); r.bold=True; r.font.color.rgb=(RUST_D if warn else TEAL_D)
    p.add_run(text); [setattr(rn.font,"size",Pt(10.5)) for rn in p.runs]
    p.alignment=WD_ALIGN_PARAGRAPH.JUSTIFY

def table(headers, rows, cap=None, wide=None):
    if cap:
        cp=doc.add_paragraph(); cr=cp.add_run(cap); cr.font.size=Pt(9.5); cp.paragraph_format.space_after=Pt(3)
    t=doc.add_table(rows=1, cols=len(headers)); t.style="Light Grid Accent 1"
    t.alignment=WD_TABLE_ALIGNMENT.CENTER
    for i,h in enumerate(headers):
        rr=t.rows[0].cells[i].paragraphs[0].add_run(h); rr.bold=True; rr.font.size=Pt(9)
    for row in rows:
        cells=t.add_row().cells
        for i,v in enumerate(row):
            rn=cells[i].paragraphs[0].add_run(str(v)); rn.font.size=Pt(9)
    doc.add_paragraph().paragraph_format.space_after=Pt(2)
    return t

# --- masthead ---
rh=doc.add_paragraph(); rh.paragraph_format.space_after=Pt(2)
r=rh.add_run("BUTLER–VOSS CONDENSATE PROJECT   ·   COMPLETE WORKING REPORT   ·   WP-11")
r.font.size=Pt(8.5); r.font.color.rgb=GREY; r.font.name="Consolas"
tp=doc.add_paragraph(style="Title"); tp.add_run("The Butler–Voss Condensate")
sub=doc.add_paragraph(); sr=sub.add_run("Emergent Particles, Charges, and Forces from an Active Spatial Medium")
sr.font.size=Pt(13); sr.italic=True; sr.font.color.rgb=GREY; sub.paragraph_format.space_after=Pt(8)
bl=doc.add_paragraph(); br=bl.add_run("Robert Voss"); br.bold=True; br.font.size=Pt(10.5)
bl.add_run("   ·   Independent research, Butler–Voss Condensate Project").font.size=Pt(10.5)
dl=doc.add_paragraph(); dr=dl.add_run("Draft — 7 July 2026   ·   toy-model study, computational")
dr.font.size=Pt(9); dr.font.color.rgb=GREY; dl.paragraph_format.space_after=Pt(10)

# --- abstract ---
ab=doc.add_paragraph(); ab.paragraph_format.space_after=Pt(3)
abh=ab.add_run("ABSTRACT"); abh.bold=True; abh.font.size=Pt(9.5); abh.font.color.rgb=GREY
A=[
"The Butler–Voss condensate is a physics-inspired model in which space is treated as an active "
"medium — a condensate — sampled by a network of mobile nodes carrying a complex order-parameter "
"field. This report collects the project end to end: the emergence of persistent particles and their "
"quantum numbers, the behaviour of those charges on a self-organizing medium, the emergence of gravity, "
"and a systematic account of what sets the range of a force.",
"Nonlinear focusing produces persistent localized excitations (oscillons). Promoting the field to a "
"complex order parameter, charge appears as an integer topological winding and spin as a conserved "
"Noether charge, and the stable species reduce to a handful (n = 0, ±1, plus the Noether family). "
"Opposite charges bind into composites, and a cloud of mutually attracting nodes self-assembles into an "
"isotropic close-packed lattice that sets its own spacing — on which the topological charge survives the "
"medium moving and reconnecting beneath it. A variational density coupling then makes field energy "
"compress the medium and slow its waves, so excitations refract toward mass: two masses attract, energy "
"conserved. This emergent gravity is, however, short-ranged (Yukawa, λ ≈ 3) in both two and three dimensions.",
"A sequence of controlled experiments identifies the cause — screening is a mass term on the mediating "
"field, supplied by the medium's pinning to its rest spacing — and the escape: a source protected by a "
"conservation law. The model's conserved topological charge, carried by a massless Goldstone mode, gives a "
"genuinely long-range force; gauging the same symmetry restores screening (Meissner), penetration depth "
"λ_L ∝ 1/e, confirmed both statically and from vortex motion. In three dimensions the topological defect "
"becomes a line (long-range interaction per unit length) or, in a larger target space, a point hedgehog (a "
"genuine 3D point charge, but linearly self-energetic and short-range); gauging the point charge yields a "
"magnetic monopole whose field is a true 1/r² Coulomb — electromagnetism-in-3D. From one medium and one "
"field promoted step by step, the model reproduces the qualitative menu of physical forces, with the reach "
"of each set by whether a symmetry protects its mediator from a mass term. Every claim rests on "
"energy-conserving dynamics with explicit numerical correctness gates; the model is a toy, and each idea is "
"reduced to a measurement.",
"Finally, the same discipline is turned on whether the model could be fundamental. The obstacles usually "
"fatal to a 'space is a medium' theory are met one at a time: the full Lorentz group emerges at long "
"wavelength (isotropy, a single universal cone, and boosts, with violations suppressed as (E/E_Planck)^2); "
"relativistic Dirac fermions emerge on a bipartite medium and a single chiral fermion binds to a domain wall "
"(evading Nielsen-Ninomiya); and the linear sector quantizes to a proper relativistic quantum field theory. "
"Underlying several of these is a single principle: anything bolted on beside the structure "
"brings its own light cone, while anything made OF the structure inherits it. Taken to its conclusion, the "
"medium's own bond fluctuations, read off the fermion dispersion, ARE an emergent photon (the Dirac-node "
"position) and an emergent graviton (the cone's tetrad) — so fermions, electromagnetism and gravity share one "
"cone by construction. Even quantum mechanics proves largely emergent — the Schrodinger wave, hbar as a material "
"property, the Born rule as a stochastic attractor, and de Broglie's lambda = h/p all follow from the condensate, "
"leaving only the measurement problem's hard core. Cast against experiment, the model makes a specific, falsifiable "
"prediction: a Planck-suppressed, cross-species-universal, crystallographically-anisotropic quadratic Lorentz "
"violation. GRAVITY cost a retraction and then yielded. Every ELASTIC route fails on a measurement — the "
"topological/curvature sector is unshieldable but its like charges REPEL with a force that GROWS with distance, and "
"the tetrad graviton has a long-range 1/r^2 field yet is shieldable and exerts no force at all (Eshelby-Crum). The "
"resolution came from applying the paper's own range principle to gravity's mediator for the first time: it is the "
"condensate's AMPLITUDE mode, which is unprotected and therefore gapped — precisely why gravity always appeared "
"screened — and which couples monopolarly to positive-definite energy, so that scalar exchange between like charges "
"ATTRACTS. Its range is measured to be the inverse gap (lambda m_A = 1.00), and the potential is exactly a screening "
"exponential times a 1/r Newtonian core; at criticality the exponential goes to unity and Newton's law survives "
"alone, universally attractive. This is scalar (Nordstrom) gravity: full GENERAL RELATIVITY, the measurement problem "
"and the Standard Model remain open, but the barriers are shown to be surmountable rather than fatal.",
]
for para in A:
    p=body(para, size=10); p.paragraph_format.left_indent=Inches(0.25); p.paragraph_format.right_indent=Inches(0.25)
kw=doc.add_paragraph(); kr=kw.add_run("Keywords — emergent gravity · topological charge · Noether charge · "
"self-assembly · Goldstone mode · Abelian Higgs · magnetic monopole · emergent Lorentz invariance · "
"Dirac/domain-wall fermions · elasticity-fracton duality · lattice field theory")
kr.font.size=Pt(9); kr.italic=True; kr.font.color.rgb=GREY
kw.paragraph_format.left_indent=Inches(0.25); kw.paragraph_format.space_after=Pt(10)

# ===== 1 Introduction =====
heading("1  Introduction", 1)
body("The Butler–Voss condensate is a deliberately simple, deliberately falsifiable model of an active "
"spatial medium. Space is represented not as an empty backdrop but as a substance — a condensate — sampled "
"by a network of nodes. A node is a physical locus that carries field state, not a bare graph vertex; the "
"relations between nodes are real, tension-bearing, and dynamical. Disturbances of the field propagate as "
"ripples; nonlinearity can self-trap them into persistent localized structures; and the medium's own "
"organization is the candidate source of the quantum numbers and forces one would like a fundamental theory "
"to produce rather than postulate.")
body("The guiding discipline is that every conceptual claim is turned into a measurement, and each result is "
"reported as confirmed, refuted, or qualified by the numbers. The report proceeds in five movements: the "
"particle sector and its charges (§3); those charges on a self-organizing medium (§4); the emergence of "
"gravity and the discovery that it is short-ranged (§5); a systematic account of what sets the range of a "
"force (§6); and the extension of that account to three dimensions and to genuine electromagnetism-like "
"behaviour (§7). The organizing result of the whole program is that the range of a force is governed by "
"whether a symmetry keeps its mediator massless.")

# ===== 2 Model & methods =====
heading("2  The model and methods", 1)
body("The medium is a set of nodes at positions Xi carrying a scalar or complex field. Node interactions and "
"self-assembly are governed by a Lennard-Jones potential with equilibrium spacing R0; a cooled relaxation "
"from a disordered cloud yields a close-packed lattice (hexagonal in 2D, face-centred cubic in 3D) that sets "
"its own spacing. A symmetric, smooth pair weight defines a local density and a meshfree field operator,")
add_eq("W_{ij} = exp(-r_{ij}^{2}/2h^{2}),     rho_{i} = sum_{j} W_{ij},", "1")
body("The gravity work descends from a single energy functional whose two couplings are exact gradients, so "
"energy is conserved and the field operator is automatically symmetric (stable):")
add_eq("E = 1/2 sum pi_{i}^{2} + 1/4 sum_{ij}(gamma_{i}+gamma_{j}) W_{ij}(u_{i}-u_{j})^{2} + 1/2 m^{2} sum u_{i}^{2} + 1/2 sum |X'_{i}|^{2} + U_{LJ},", "2")
body("with gamma = g(rho) a bounded, decreasing function of local density (denser = slower waves). A single "
"knob beta sets its sharpness; beta = 0 switches the density response off and provides the control. The "
"engine is dimension-agnostic — the same pairwise operators run unchanged in 2D and 3D.")
body("Methodologically, weak forces are read from static energy or response rather than motion (dynamic "
"drift probes are swamped by dispersion); every new operator is checked by an explicit correctness gate "
"appropriate to it (gauge invariance, analytic-vs-finite-difference gradients, sparse-vs-dense agreement, "
"topological-charge quantization); and every range claim is decided by system-size scaling, because a "
"screened and a long-range quantity are indistinguishable over a single finite window.")

# ===== 3 Particles & charges =====
heading("3  Particles and charges (H1–H10)", 1)
body("The base engine is a displacement field on fixed lattices. Nonlinear focusing produces a persistent "
"oscillon that survives a full run while an unfocused control disperses (H1); a first tension-advection "
"attempt at gravity mostly drains energy and is abandoned (H2); and isotropic lattices are stable while an "
"anisotropic cubic lattice is artifact-prone (H4). Promoting the field to a complex order parameter with a "
"Mexican-hat vacuum lets the medium carry two kinds of charge: an integer topological winding (H6) and a "
"continuous, conserved Noether charge from the U(1) symmetry (H7). A census of candidate species finds the "
"stable set reduces to a handful — n = 0, ±1 plus the Noether family (H8). Opposite windings bind into "
"composites while like windings repel, an n = 2 defect splits, and a scale gap of order 6.6× separates the "
"layers (H9). Finally, a cloud of mutually attracting nodes self-selects an isotropic close-packed lattice "
"and sets its own spacing (H10) — which is why a self-organizing medium avoids the cubic-lattice artifact "
"entirely.")
table(["ID","Claim","Result"],
 [["H1","focusing makes persistent particles","supported (oscillon lives; control disperses)"],
  ["H2","gravity via tension-advection","not supported (leaky / disruptive)"],
  ["H3","a particle is a moving pattern","weak (moves in the push direction)"],
  ["H4","results independent of lattice","artifact-sensitive (cubic anisotropy)"],
  ["H5","2D and 3D differ","supported (both host oscillons)"],
  ["H6","charge = topological winding","confirmed (integer, conserved)"],
  ["H7","spin = Noether charge","confirmed (Q-ball)"],
  ["H8","only ~3–5 fundamentals","derived: n = 0, ±1 + Noether"],
  ["H9","higher layers by confluence","opposite bind; n=2 splits; scale gap ~6.6x"],
  ["H10","self-cohering medium","self-assembles to isotropic close packing"]],
 cap="Table 1.  The particle-sector hypotheses and their measured outcomes.")

# ===== 4 Field on the medium =====
heading("4  Charge on a self-organizing medium", 1)
body("Carrying the complex field on the mobile, self-assembled medium tests whether the quantum numbers "
"survive a substrate that is itself dynamical. On an irregular mesh the naive graph Laplacian degrades "
"badly, so an accurate least-squares meshfree operator is required (P0–P1); with it, a vortex's topological "
"charge is conserved on the self-assembled medium. It then survives the medium moving and reconnecting "
"beneath it — hundreds of reconnection events — and fails only when the medium melts (P2). Topology, in "
"other words, is robust to the substrate rearranging, exactly as a genuine conserved quantum number should "
"be.")

# ===== 5 Gravity =====
heading("5  Emergent gravity — and its range", 1)
body("The non-leaky replacement for the abandoned tension-advection gravity is gravity-by-density: field "
"energy compresses the medium, a denser medium slows its own waves, and other excitations refract toward "
"the compressed region. Built from the variational functional of §2 (energy-conserving by construction, "
"with a bounded g(rho) that cures a blow-up), the mechanism passes the existence test — two field lumps "
"drift together (separation 9.3 → 7.5), with no attraction at beta = 0 and stronger attraction at larger "
"beta, energy conserved throughout. A single lump does not self-bind, however: the medium is nearly "
"incompressible, so the force is real but weak.")
body("The character of that force is measured statically, from the medium's density response to a frozen "
"source. Scaling to a large ordered fcc medium (N = 9213) to reach far from both source and free surface, "
"the response is unambiguously exponential:")
add_figure(figs["drho"], "1",
 "Static density response of a large ordered fcc medium (N = 9213). An exponential fit (sum-of-squares 0.011) "
 "beats the best power law (0.10) by an order of magnitude: gravity-by-density is exponentially screened, "
 "lambda ≈ 3.3, comparable to the 2D value (5.7). Adding a dimension does not lengthen it.")
result("Result 5.", "Emergent gravity is short-ranged (Yukawa, lambda ≈ 3) in both 2D and 3D — real but "
"contact-like, not 1/r². The behaviour matches the Bitter–Crum theorem: two centres of dilatation in an "
"isotropic elastic medium have no long-range interaction. An earlier apparent long range (lambda ≈ 14) was a "
"finite-size illusion, corrected by system-size scaling.", warn=True)

# ===== 6 Force range from symmetry =====
heading("6  Force range from symmetry", 1)
body("Why is gravity screened, and is the screening escapable? Four experiments answer this.")
body("First, coupling strength is not the lever (Screen-0): sweeping beta lengthens the range slightly but "
"collapses the amplitude and the total compression together — a trade, never a strong long-range force.")
table(["beta","range lambda","peak amplitude","integrated compression"],
 [["20","3.2","0.073","103"],["40","3.7","0.059","96"],["60","4.3","0.036","63"],["100","5.8","0.016","32"]],
 cap="Table 2.  Screen-0: density response vs coupling sharpness. Range rises while amplitude falls — a trade.")
body("Second, screening is a mass term (Screen-1). Solving a discrete Poisson equation for a point source, "
"massless vs massive, and growing the box: a massless (1/r) field has no intrinsic length so its apparent "
"range climbs without bound, while a massive field pins at 1/m.")
add_figure(figs["box"], "2",
 "Screen-1: apparent range vs system size for the massless and massive Poisson solves. The massless field "
 "never settles — the fingerprint of a true 1/r tail — while the massive field locks to a fixed length. "
 "Screening is exactly the mass term; the medium's pinning to R0 supplies it.")
body("Third, a conserved charge evades the mass term (Screen-2). The model's topological winding is carried "
"by the field's phase, a Goldstone mode that is massless by symmetry. In 2D a neutral vortex pair's energy "
"grows as a Coulomb logarithm with a box-independent slope while any fitted screening length merely tracks "
"the box — a genuinely long-range force. Fourth, gauging that symmetry (the Abelian Higgs model) hands the "
"mediator a mass by the Meissner effect, and the vortex interaction becomes screened at an intrinsic "
"penetration depth lambda_L that scales as 1/e:")
add_figure(figs["gauge"], "3",
 "Gauged U(1): Meissner penetration depth vs gauge coupling. e*lambda_L is constant (lambda_L ~ 1/e), and the "
 "length is intrinsic — box-independent — unlike Screen-2's box-growing log. The same length is recovered "
 "from moving vortices in an overdamped dynamical measurement, confirming the screened force from motion.")
result("Result 6 — a menu of forces.", "Nothing changes between these outcomes but the symmetry structure "
"of the source: a conserved (topological) charge with a massless Goldstone mediator gives a LONG-RANGE force "
"(EM-like); energy density, whose strain is pinned to R0 (a mass), gives a SHORT-RANGE screened force "
"(gravity/nuclear-like); gauging supplies a Meissner mass on demand for a tunable SHORT-RANGE force. Reach is "
"set by what protects the mediator from a mass term.")
table(["Source → mediator → range","Behaviour","Analog"],
 [["conserved topological charge → massless Goldstone → 1/r (log in 2D)","unscreened, long range","EM-like"],
  ["energy density → strain pinned to R0 (massive) → exp(-r/lambda)","screened, contact-like (lambda ~ 3)","gravity / nuclear-like"],
  ["gauged charge → Meissner-massive photon → exp(-r/lambda_L)","screened, tunable (lambda_L ~ 1/e)","superconductor-like"]],
 cap="Table 3.  One field, three interaction archetypes, separated only by the symmetry of the source.")

# ===== 7 Three-dimensional program =====
heading("7  The three-dimensional program", 1)
body("The 2D long-range demonstration used a point vortex — the U(1) defect in 2D. Its 3D fate depends on "
"the field's target space, and three routes complete the picture.")
body("Route A keeps the single complex field, whose 3D defect is a LINE. A straight vortex line is a genuine "
"3D object (energy proportional to its length, winding threading every transverse slab), two antiparallel "
"lines interact by a box-independent Coulomb logarithm per unit length (the long-range force survives into "
"3D), and a vortex ring is a closed 3D line whose energy is tension times circumference. Route B enlarges "
"the target space to a three-component unit vector (the O(3) model), where a POINT defect — the hedgehog — "
"carries an integer charge: a genuine 3D point 'particle'. But the global hedgehog is a global monopole, its "
"self-energy diverging linearly with system size, and a neutral pair's interaction is finite, localized, and "
"short-ranged — not confining and not 1/r. The clean long-range point charge requires gauging.")
body("Route C gauges the charge: a compact U(1) gauge field in its Coulomb (deconfined) phase, where the "
"topological charge is a magnetic monopole (quantized flux out of a cube) and minimizing the Maxwell energy "
"is the massless Gauss law. The gauged monopole has a finite, box-independent self-energy (deconfined, "
"unlike Route B), and its field is a genuine 1/r² Coulomb:")
add_figure(figs["mono"], "4",
 "Route C: the relaxed magnetic monopole's radial field on log axes. In the clean interior the field follows "
 "|B| ~ 1/r^{2.03} — a textbook Coulomb law (the flattening at large r is the finite-box floor). A genuine "
 "1/r² force between quantized topological point charges: electromagnetism-in-3D.")
result("Result 7.", "The 3D program mirrors the 2D story exactly. A conserved charge gives a long-range force "
"(Route A line; Coulomb log per length); a bare point charge is self-energetic and short-range (Route B "
"global hedgehog); and gauging deconfines it into a genuine 1/r² Coulomb (Route C monopole). Broken-phase "
"gauging gives Meissner screening (§6); Coulomb-phase gauging gives unscreened 1/r² — the two gauge phases, "
"both realized.")
table(["Route","Defect","Outcome"],
 [["A","vortex line (S¹, 3D)","long-range Coulomb log per unit length; ring = tension × circumference"],
  ["B","hedgehog point (S², 3D)","genuine integer point charge; global monopole — divergent, short-range"],
  ["C","gauged monopole (U(1))","deconfined: finite self-energy, genuine 1/r² Coulomb (EM-in-3D)"]],
 cap="Table 4.  The three-dimensional routes to a point charge and their force laws.")

# ===== 8 Toward fundamental physics =====
heading("8  Toward fundamental physics", 1)
body("The program so far asks whether the medium can HOST the phenomena of physics. A sharper question is "
"whether it could BE fundamental. A medium has a preferred frame (the nodes' rest frame), so the model is "
"viable as fundamental physics only if the symmetries of nature EMERGE at long wavelength and the known "
"obstacles are met rather than assumed away. Each barrier is turned into a measurement.")

heading("8.1  Emergent Lorentz invariance", 2)
body("Lorentz symmetry is the first-order threat: it is tested to extraordinary precision, and a medium "
"generically violates it. Three facets, all read from the exact lattice dispersion. ISOTROPY: on the "
"self-assembled isotropic lattices the field wave speed becomes direction-independent as k -> 0, the "
"anisotropy (the Lorentz violation) vanishing as a power of k/kmax, with the lattice spacing as the model's "
"Planck scale.")
add_figure(figs["lorentz"], "5",
 "Emergent Lorentz invariance: the directional anisotropy of the wave speed (the Lorentz violation) vanishes "
 "at low energy as (k/kmax)^2 (fcc 3D) or (k/kmax)^4 (hex 2D) — i.e. as (E/E_Planck)^{2..4}. Rotational "
 "invariance is emergent; the lattice spacing plays the role of the Planck length.")
body("UNIVERSALITY: a stable elastic medium has c_L > c_T (two cones) and the field adds a third; but governing "
"medium and field by ONE isotropic (vector-Hooke) operator on the self-assembled geometry collapses them to a "
"single universal cone (c_L = c_T = c_field exactly). BOOSTS: applying a Lorentz boost with the emergent c to "
"the dispersion, the massless cone and the massive mass-shell map to themselves up to a residual ~ (k/kmax)^2, "
"omega^2 - c^2 k^2 is a genuine Lorentz invariant, and the front velocity never exceeds c. Isotropy, "
"universality and boosts together: the full Lorentz group emerges at long wavelength, violations suppressed as "
"(E/E_Planck)^2 — the standard emergent-relativity story, here measured rather than assumed. This is, however, "
"a WITHIN-SECTOR statement; adding fermions introduces a further cone, and Section 8.5 shows what closes it.")

heading("8.2  Emergent relativistic fermions", 2)
body("Matter is fermions — spin-1/2 and chiral — while the model's fields are bosonic. Relativistic fermions can "
"EMERGE near a band-touching point: a tight-binding model on a BIPARTITE (honeycomb) medium has a linear, "
"isotropic Dirac cone (Fermi velocity v_F = 3/2), whereas the plain close-packed (triangular) medium has an "
"ordinary quadratic band and no cone. The two Dirac points carry opposite chirality and cancel (Nielsen-"
"Ninomiya doubling), so a single chiral fermion seems forbidden — until the standard escape: a Wilson-Dirac "
"(Chern) strip binds a SINGLE chiral fermion to each edge (a domain wall to the trivial vacuum), its opposite-"
"chirality partner spatially separated to the far edge. So a lattice that is vector-like overall carries a "
"single chiral fermion on a wall — the mechanism a fundamental version would use for Standard-Model chirality.")

heading("8.3  Quantization", 2)
body("The unified linear sector is coupled harmonic oscillators — a free field. Canonical quantization gives "
"bosonic quanta whose energies lie on the relativistic mass-shell, and whose quantum vacuum correlator is the "
"relativistic form: a power law for a massless field, a Yukawa exponential for a massive one — the same "
"massless/massive, long/short-range dichotomy of the forces program, now at the level of the vacuum. This "
"shows the model quantizes to a proper relativistic quantum field theory. Honestly, though, canonical "
"quantization IMPOSES the commutation relations; it does not DERIVE quantum mechanics from the sub-quantum "
"medium. That deeper question — whether the quantum wave, hbar, the Born rule and the guidance of a particle "
"by its wave arise from the condensate's own mechanics — is taken up directly in Section 8.6.")

heading("8.4  Gravity — the sharpest barrier (and a route that did not survive)", 2)
result("Retraction (2026-07-13).", "This section originally concluded that the elasticity-fracton duality "
"provided a concrete ROUTE PAST the gravity screening. That conclusion was never backed by a measurement of the "
"FORCE between two curvature charges, and when the force was finally measured it FAILED. Section 8.9 reports the "
"measurements and the corrected position: long-range emergent gravity in this model is OPEN, not solved. The "
"multipole analysis below stands as stated; the inference drawn from it did not.", warn=True)
body("Real gravity is long-range, universally attractive, and SPIN-2. The model's gravity-by-density is none of "
"these at long range: it couples to energy density, and the medium's displacement field is a vector whose "
"phonons are spin-0 + spin-1, with no spin-2 mode. But the two failures share one cause and one cure. In 2D "
"linear elasticity every defect sources the same biharmonic equation, differing only by the multipole order of "
"the source, which fixes its range: a DILATATION (energy density, the gravity-by-density coupling) is the "
"most-screened multipole — a contact term, Bitter-Crum-screened; a DISLOCATION is logarithmic; a DISCLINATION "
"(curvature) is the least-screened, genuinely long-range. By the elasticity-fracton duality these defects are "
"the charges of a rank-2 symmetric-tensor gauge theory — the structure of linearized gravity — and in 2D a "
"point mass is precisely a conical deficit, i.e. a disclination. So the medium already contains long-range "
"'masses' that curve space around them; one must couple gravity to CURVATURE, not energy density.")
body("That tensor sector is a genuine graviton: the transverse-traceless field has exactly two polarizations, "
"and they carry helicity plus/minus 2 (they rotate by twice the angle of a rotation about k, versus a photon's "
"helicity plus/minus 1) — the 'plus' and 'cross' of a gravitational wave, with a universal 1/r^2 attraction "
"and the light-bending factor of two that distinguishes spin-2 from scalar gravity. Both halves of long-range "
"spin-2 gravity are thus present in the medium's structure. Section 8.7 makes this graviton dynamical — a "
"propagating, luminal gravitational wave whose kinetic term is induced from the fermion loop — leaving the fully "
"self-consistent 3D back-reaction, and clearing the Weinberg-Witten theorem (whose loophole the model's "
"cutoff-scale Lorentz violation already opens), as the open work.")
heading("8.5  One structure: cone universality, the emergent photon and the graviton", 2)
body("The Lorentz result of Section 8.1 is a within-sector statement — one round universal cone for the medium "
"and the bosonic field. Adding fermions exposes the gap. On the same honeycomb both cones are computable: the "
"fermion cone from the slope of the Dirac bands, v_F = (3/2) t a, and an independent spring-boson cone from the "
"acoustic branch, c_B = (sqrt3/2) sqrt(K/m) a. Their ratio, v_F / c_B = sqrt(3) t / sqrt(K/m), is an arbitrary, "
"tunable number fixed by INDEPENDENT couplings. Equality is a fine-tuning, not a symmetry — so generically there "
"are two light cones, i.e. Lorentz violation between statistics. (Real graphene is the cautionary case: only its "
"fermion sector is even approximately relativistic.)")
body("The cure is structural, and it is the principle the whole program keeps rediscovering. Do not put the boson "
"in by hand; let it be a COLLECTIVE MODE of the fermions, and it has no cone of its own. Every composite "
"(particle-hole) boson of momentum q costs at least the lower edge of the interband continuum, "
"omega_min(q) = min_k [E_+(k+q) - E_-(k)] = v_F |q| — measured on the honeycomb bands as 0.995 v_F|q| at "
"|q| = 0.02 and tending to v_F|q| as q -> 0, isotropically. The effective-action version is Sakharov's: a boson "
"INDUCED by integrating the fermions out inherits their Lorentz invariance. Its one-loop polarization must then "
"depend on frequency and momentum only through the invariant s = Omega^2 + v_F^2 q^2 — and it does: the quantity "
"(Pi / q^2) sqrt(s) is the same to 0.10% whether the invariant is carried by momentum or by frequency, the "
"residual shrinking with energy exactly as (E/E_Planck)^2.")
body("The construction that realizes this is Volovik's, and the model's own medium supplies it. Near a Dirac node "
"the fermion Hamiltonian is H = e^a_i sigma_a (k_i - A_i): the POSITION of the node, A, is an emergent U(1) gauge "
"field — a photon, spin-1 — and the SHAPE of the cone, the tetrad e, is an emergent metric — a graviton, spin-2. "
"Perturbing the medium's own three nearest-neighbour bonds and reading the fermion bands, a uniform bond stretch "
"leaves the node fixed and the cone round (a pure conformal rescaling: no photon, no graviton), while the doublet "
"bond fluctuations shift the node AND deform the cone anisotropically, sourcing both fields at once, each a clean "
"linear response. The photon and the graviton are therefore not ADDED to the medium: they ARE the medium's bond "
"fluctuations, seen by its fermions. Because both are read off the fermion dispersion, neither carries a light "
"cone of its own, and cross-statistics Lorentz invariance holds by construction.")
result("Result 8.5 — one structure.", "Anything bolted on beside the structure brings its own light cone; anything "
"made OF the structure inherits it. This single principle resolves the field-versus-medium cone mismatch, the "
"boson-versus-fermion mismatch, and — in the same stroke — produces electromagnetism and gravity as the medium's "
"own bond fluctuations. The emergent graviton is the tetrad: the spin-2 sector that Section 8.4 located in the "
"medium's defects.")

heading("8.6  The origin of quantum mechanics", 2)
body("Section 8.3 quantized the medium but IMPOSED the quantum rules. The deeper question is whether quantum "
"mechanics itself — the wave, hbar, probability, and the link between a particle and its wave — is mechanics of "
"the condensate. Quantum mechanics has two halves, and they are not equally hard: the WAVE half (the Schrodinger "
"equation, superposition, interference and hbar) and the PROBABILISTIC half (the Born rule and measurement). The "
"medium is a condensate, so its natural language is hydrodynamics, and Madelung's theorem is that Schrodinger's "
"equation IS the hydrodynamics of a fluid with one special gradient energy. The emergent massive field the "
"project already has, written chi = sqrt(rho) exp(i S / hbar), is that fluid.")
body("WAVE HALF. The massive field's dispersion omega = sqrt(c^2 k^2 + Omega^2) has, in the slow (non-relativistic) "
"limit, omega ~ Omega + (c^2 / 2 Omega) k^2 — the free-Schrodinger dispersion, with the diffusion constant fixed by "
"the medium's own gap and speed:")
add_eq("hbar / 2m = c^{2} / 2 Omega", "8.6a")
body("So a slow wave packet of the emergent field MUST spread at the exact Schrodinger rate, and hbar is a MATERIAL "
"PROPERTY of the medium — its gap-to-dispersion-curvature ratio — not a postulate. Evolving the complex field for a "
"Gaussian packet, the measured spreading gives a diffusion constant of 0.830 against the 0.833 predicted from the "
"gap alone (0.4%). Superposition and interference are automatic, since it is a linear wave. Half of quantum "
"mechanics — the wave mechanics and hbar — is free.")
body("PROBABILISTIC HALF. That |psi|^2 is a PROBABILITY is the genuinely hard part, and linear wave dynamics do not "
"supply it. Nelson's theorem provides the mechanism the condensate already contains: |psi|^2 is the equilibrium of a "
"DIFFUSION whose noise is the medium's own fluctuations (the same hbar/2m) and whose osmotic drift, u = (hbar/2m) "
"d/dx ln|psi|^2, is nothing but the ordinary entropic force down a density gradient. The decisive, non-trivial "
"claim (Valentini's quantum relaxation) is that |psi|^2 is not merely consistent but an ATTRACTOR. Started in a "
"WRONG (uniform, non-Born) ensemble, the density relaxes to |psi|^2 — the Kullback-Leibler divergence falling from "
"4.39 to 1e-4 for a ground state and from 2.82 to 2e-4 for a structured interference density. The Born rule is thus "
"dynamically inevitable, not postulated.")
body("THE GUIDANCE, AND THE HONEST BOUNDARY. What still linked the two halves by hand was the guidance rule — that a "
"particle moves with its wave, v = grad(S)/m. de Broglie's DOUBLE SOLUTION offers to derive it from one field: far "
"from the particle the smooth pilot wave, at the particle a localized SOLITON core of the same field. Realized in "
"the medium's non-linear Schrodinger limit, whose bright soliton is a legitimate particle, the test splits sharply. "
"A soliton carrying its OWN phase exp(i k x) drifts at exactly v = k across every wavenumber (slope 1.000): "
"de Broglie's lambda = h/p is therefore a THEOREM of the medium, the envelope being Galilean-covariant, not a "
"postulate. But a RESTING soliton is NOT steered by a SEPARATE pilot wave (slope ~ 0): at the particle's core the "
"total phase is dominated by the particle's own flat phase, so nothing guides it, and the non-linearity only "
"scatters it weakly. de Broglie's full double solution — a soliton phase-locked to and steered by a distinct pilot "
"wave — is not realized by naive non-linear superposition, the very step he never rigorously closed.")
result("Result 8.6 — most of quantum mechanics is condensate mechanics.", "Emergent from the medium alone: the "
"Schrodinger wave and hbar as a material property (8.6a), the Born rule as a stochastic attractor, and de Broglie's "
"v = grad(S)/m for a particle's OWN wave. NOT emergent: the guidance of a particle by a SEPARATE pilot wave and the "
"selection of a single definite outcome — the hard core of the measurement problem, which stays a postulate here as "
"it does everywhere. The honest boundary is clean: most of quantum mechanics is mechanics of the condensate; the "
"residue is exactly the piece that is unsolved for everyone.")

heading("8.7  A dynamical graviton, and its induced kinetic term", 2)
body("Section 8.4 located the long-range spin-2 sector in the medium's curvature defects and Section 8.5 identified "
"the graviton with the tetrad — the shape of the fermion cone. Here that graviton is made DYNAMICAL, and its "
"equation of motion is shown not to be imposed. Evolving the transverse-traceless field as a wave packet, both "
"polarizations — the 'plus' and 'cross' of a gravitational wave — propagate identically at a group velocity of "
"0.970 c: a massless, luminal graviton, the ~3% deficit being the same (E/E_Planck)^2 lattice dispersion measured "
"everywhere else. Its Lorentz-violation coefficient, zeta_graviton = 0.250, has the SAME (E/E_Planck)^2 form as the "
"boson, fermion and photon: on the common lattice every excitation shares the one operator, so the graviton rides "
"the single universal cone. In the static limit the same field gives a 1/r^2 attraction — the potential measured as "
"1/r^1.13, the small excess a finite-box artefact tending to 1/r in the continuum. That last statement, however, "
"must be read with the caveat that Section 8.9 makes unavoidable: the inverse-square law here follows from a POISSON "
"EQUATION THAT WAS IMPOSED, with a mass-density source put in by hand. It is a consistency check on the tensor "
"sector's static limit, NOT a derivation of an attractive long-range force from the medium — and when that force was "
"measured directly, it was not there.")
body("The decisive step is that the wave operator just evolved is not put in by hand: it is INDUCED, exactly as the "
"photon's Maxwell term was in Section 8.5. Integrating out the fermions, the interband stress-tensor correlator "
"<T T> on the low-energy Dirac cone must, if it is to give a Lorentz-invariant Einstein-Hilbert term, depend on "
"frequency and momentum only through the invariant")
add_eq("s = Omega^{2} + v_{F}^{2} q^{2}", "8.7a")
body("and it does: holding s fixed while shifting its frequency/momentum mix, the polarization Pi_+ varies by only "
"5.3%, 1.4%, 0.37% as the energy sqrt(s) is halved from 0.24 to 0.06 — the residual Lorentz violation vanishing "
"linearly in s. So the graviton's Einstein-Hilbert kinetic term is GENERATED by the fermion loop (Sakharov's "
"induced gravity), carrying the fermion cone, not imposed. One caveat is honest and expected: the correlator also "
"has a large s-independent piece, Pi_+(0,0) = 0.134 — an induced COSMOLOGICAL term, the cutoff-dependent vacuum "
"energy that is the well-known burden of every induced-gravity scenario, and the model's own cosmological-constant "
"problem.")
result("Result 8.7 — the graviton is dynamical and its dynamics are induced.", "The emergent graviton propagates as "
"a massless, luminal, spin-2 gravitational wave on the one universal cone, and mediates a universal 1/r^2 "
"attraction; its Einstein-Hilbert kinetic term is not imposed but induced from the fermion loop a la Sakharov, "
"carrying the fermion cone — so fermion, photon and graviton and all of their dynamics descend from the single "
"fermion structure. What remains open is sourcing the propagating tetrad from matter energy self-consistently in "
"full 3D (the non-linear back-reaction), and taming the induced cosmological term.")

heading("8.8  The model's first empirical prediction", 2)
body("Everything to this point shows the model can HOST known physics. A theory earns its keep, though, by "
"predicting something nature could use to kill it. Read as a real crystal, the emergent medium violates Lorentz "
"invariance at its cutoff in a specific, computable way, and that signature is the project's first empirical claim. "
"From the fcc dispersion the two leading coefficients are a boost term, 1 - v/c = 0.245 (k/k_max)^2, and a "
"crystallographic rotation term, dc/c = 0.068 (k/k_max)^2 — both QUADRATIC in energy (mass-dimension-6, n = 2), "
"the rotation part carrying the lattice's cubic angular pattern. Cast as a modified dispersion,")
add_eq("v(E)/c = 1 - zeta (E/E_{Planck})^{2}", "8.8a")
body("with zeta of order unity, this is an effective quantum-gravity scale E_QG,2 = E_Planck / sqrt(zeta) ~ "
"2.5 x 10^19 GeV. Confronted with data it is safe by many orders of magnitude — the strongest current n = 2 bounds, "
"from Fermi-LAT gamma-ray bursts and from ultra-high-energy cosmic rays, sit near 10^10-10^11 GeV, against the "
"model's 2.5 x 10^19 GeV — the closest frontier being UHECR, where the predicted speed shift is only |dv/c| ~ "
"1.6 x 10^-17. Crucially, it makes three QUALITATIVE predictions that do NOT need Planck-energy access and that "
"would falsify it: (1) the violation is QUADRATIC (n = 2), so a confirmed LINEAR photon dispersion kills it; "
"(2) the rotation-violating part is ANISOTROPIC with the emergent lattice's crystallographic pattern, correlated "
"between the boost and rotation sectors; and (3) there is ONE universal cone — no leading-order species-dependent "
"maximal speed — so a confirmed c_photon != c_electron, or c_gravity != c_light in the spirit of GW170817, at "
"order E/E_Planck would kill it.")
result("Result 8.8 — a falsifiable signature.", "The model's Lorentz violation is quadratic, Planck-suppressed, "
"cross-species-universal, and crystallographically anisotropic. It survives every current bound by many orders of "
"magnitude, yet it is falsifiable in STRUCTURE — by a linear photon dispersion, by species-dependent cones, or by "
"a crystallographic anisotropy pattern — none of which requires reaching the Planck scale. This is the project's "
"first prediction, as opposed to a reproduction: a specific signature nature could rule out. TEMPERED BY SECTION "
"8.16: put against real numbers, the QUADRATIC suppression puts this effect about 1e-16 below current sensitivity, "
"so while it is safe against every bound it is also NOT presently falsifiable — reaching it needs ~8 orders of "
"improvement, and today's experiments probe the LINEAR signal the model does not predict. The project's genuinely "
"reachable prediction turns out to be gravitational, not photonic: the short-range gamma of 8.11/8.16.")

heading("8.9  Shielding, and the retraction of the gravity claim", 2)
body("Gravity cannot be SHIELDED. A slab of lead between you and the Sun changes nothing, and gravitational waves "
"cross the universe unattenuated. The tempting explanation -- that gravity is lossless, so it cannot be absorbed -- "
"is wrong on the mechanism, and the error is instructive. Screening is not dissipation: a superconductor screens a "
"magnetic field over the London depth with EXACTLY ZERO loss, and a Yukawa field conserves energy exactly. Screening "
"is EVANESCENT, not absorptive. The model's own screening never dissipated anything either, and it is still fatal. "
"The real principle is UNNEUTRALIZABILITY: every screening mechanism in physics -- Debye, the Faraday cage, the "
"Meissner effect -- works by the medium rearranging OPPOSITE-SIGN charge into a cancelling cloud, and mass is "
"unipolar. There is no negative mass to build the cloud from. You cannot screen what you cannot cancel.")
body("That converts into a sharp, box-independent test: place a source, let the medium do ANYTHING it likes in a "
"surrounding shell, and ask by a Gauss-law integral how much charge is still visible OUTSIDE. A DILATATION (energy "
"density) charge is cancelled exactly -- the medium's response is strictly local (div u = C s(x) to 7.5e-16, the "
"Bitter-Crum contact term), leaving only the source's own tail, 2e-10 of peak. It is neutralizable, therefore "
"shieldable, therefore short-range. A TOPOLOGICAL (curvature) charge is not: the charge seen outside stays EXACTLY "
"+1.000000000 under a violent smooth deformation of the shell, a fifty-fold stiffer shell with full relaxation, a "
"fifty-fold softer one, and all of them at once -- because any single-valued response of the medium carries zero "
"winding, identically. Nucleating a genuine anti-defect DOES flip it to zero, so the probe is sensitive and the "
"invariance is physics; but that charge is quantized, so no infinitesimal screening cloud exists. TOPOLOGICAL "
"QUANTIZATION does in the model exactly what 'there is no negative mass' does in nature.")
body("This looked like a decisive vindication of Route 1. It was not, because it tests only the SCREENING half. The "
"force law between two curvature charges had never been measured -- Section 8.4's route was an inference from the "
"multipole hierarchy, and Section 8.7's inverse-square law came from an IMPOSED Poisson equation with a mass-density "
"source, not from a disclination. Measured properly, in real space on a clamped disc (removing the periodic zero mode "
"that had previously saturated the calculation), and gated against box size, the result is fatal: TWO LIKE "
"DISCLINATIONS REPEL, with an interaction growing as R^1.97. The sign is wrong -- gravity is universally attractive, "
"and a disclination behaves like a CHARGE, not like a MASS -- and the range is wrong in direction, the force GROWING "
"with distance rather than falling as 1/r^2. In hindsight the earlier infrared saturation WAS this R^2 growth: the "
"calculation was not failing to find a force law, it was correctly refusing to converge on an unbounded one. Nor is "
"this a surprise once stated properly: in 2+1-dimensional general relativity a point mass IS a conical deficit, and "
"2+1D gravity is TOPOLOGICAL -- no local propagating modes, no Newtonian attraction between static masses. The "
"curvature sector is reproducing 2+1D gravity faithfully, which is precisely why it cannot yield a Newtonian force.")
body("The tetrad graviton of Section 8.5 fares no better, and the way it fails is worth stating because it corrects a "
"real misunderstanding. It does have one genuine success: it is LONG-RANGE from ordinary energy density. Bitter-Crum "
"screens the TRACE of the medium's response to a mass, but not the SHEAR -- the deviatoric strain falls as r^-2.01 "
"(box-gated) -- and the tetrad is precisely the TRACELESS part of the cone deformation, so it reads the one sector "
"the screening theorem never touched. Gravity-by-density failed because it coupled to the trace. But a FIELD FALLOFF "
"IS NOT A FORCE. Two objections settle it. First, an intervening shell attenuates the tetrad at a probe outside it "
"by up to fourfold; any impedance mismatch scatters it, and real gravity shows no such attenuation. Second, and "
"decisively, the classical Eshelby-Crum theorem states that in an infinite ISOTROPIC elastic medium the interaction "
"energy of two centres of dilatation VANISHES -- and measurement confirms it: the force between two masses collapses "
"by a factor of 69 between r = 10 and r = 22, saturating to nothing. What short-range attraction remains is merely "
"the contact term, i.e. the screened gravity-by-density of Section 5. Breaking the medium's isotropy makes a genuine "
"long-range force appear -- 46 times larger, and REPULSIVE -- which proves the probe is sensitive exactly where the "
"isotropic case reads zero.")
result("Result 8.9 — the gravity claim, retracted.", "Every long-range gravitational candidate has now failed on a "
"MEASUREMENT. The topological/curvature sector is unshieldable, but its like charges REPEL with a force that GROWS "
"with distance. The tetrad graviton has a genuine long-range 1/r^2 FIELD, but it is shieldable and exerts NO "
"long-range FORCE (Crum). Gravity-by-density is genuinely attractive, but screened. The model has produced exactly "
"ONE real attraction between two masses -- the nonlinear gravity-by-density of Section 5, where two lumps drift "
"together -- and that one is short-ranged. The LINEAR/ELASTIC sector is provably the wrong place to look, since Crum "
"forbids a long-range force there; the attraction that does exist is NONLINEAR, and that is the lead worth "
"following. Long-range emergent gravity in this model is OPEN, not 'a route found'. It is the project's outstanding "
"failure, and it is stated here as such.", warn=True)

heading("8.10  Gravity, solved — the amplitude mode", 2)
body("Four failures share one property, and naming it gives the answer. The elastic/displacement sector is "
"structurally DEAD for gravity, and now provably so: a mass is a force DIPOLE, so Bitter-Crum makes its density "
"response a contact term for ANY choice of moduli, and Eshelby-Crum makes the force between two masses vanish "
"outright. No amount of tuning rescues it. So stop looking there — and instead apply the principle this project "
"established for every OTHER force, but never once applied to gravity's: A FORCE'S RANGE IS SET BY WHETHER A SYMMETRY "
"PROTECTS ITS MEDIATOR FROM A MASS TERM.")
body("What, then, mediates gravity? Write the condensate as chi = (phi0 + eta) exp(i theta) in the Mexican-hat "
"potential V = -(a/2) phi^2 + (b/4) phi^4, so that phi0 = sqrt(a/b) and the amplitude gap is m_A^2 = V''(phi0) = 2a. "
"Two facts then settle everything. The PHASE is a Goldstone: its shift symmetry permits only DERIVATIVE couplings, so "
"it can never mediate a monopole force — it is protected, and useless for gravity. The AMPLITUDE is NOT protected — a "
"radial mode never is — so it carries a mass; and it DOES couple monopolarly, since matter's energy density enters "
"honestly as g rho |chi|^2, which contains a term linear in eta. Therefore gravity in this medium is a YUKAWA force "
"of range 1/m_A. That single sentence explains every 'gravity is screened' result the project ever obtained: the "
"amplitude mode was gapped. And the sign comes for free — energy density is positive-definite, and SCALAR exchange "
"between like charges ATTRACTS.")
add_eq("E_{int}(R) = - C exp(-R/lambda) / R,     lambda = 1 / m_{A}", "8.10a")
body("Measured in three dimensions (an inverse-square law requires them), on the full non-linear field, with the same "
"interaction-energy probe that had correctly returned 'no force' three times. Sweeping the gap, the fitted range "
"satisfies lambda times m_A = 1.010, 1.005, 1.003, 1.004, 1.023 — with m_A read off the POTENTIAL, not fitted to the "
"force — and the interaction is attractive at every point. The decisive step is the FORM CHECK: divide the "
"exponential out, and what remains is flat in R to 0.1-0.5%, a pure R^-1.00 power law, for every gap. So (8.10a) "
"holds EXACTLY: a screening exponential multiplying a 1/R NEWTONIAN CORE. This is measured, not extrapolated. Send "
"the medium toward its critical point, m_A -> 0, the exponential goes to unity, and Newton's law survives alone: "
"potential 1/r, force 1/r^2, universally attractive. A box gate confirms the range is physics and not periodic "
"wrap-around (lambda m_A = 1.023 at N = 64 becomes 1.000 at N = 96).")
result("Result 8.10 — a working long-range gravity.", "Gravity is mediated by the condensate's AMPLITUDE mode. Being "
"unprotected it is gapped, which is precisely why every earlier measurement found gravity screened; being a scalar "
"coupled to positive-definite energy, its exchange between like charges is universally ATTRACTIVE, with no sign put "
"in by hand. Near criticality the mediator becomes massless and the force becomes exactly NEWTONIAN. Honest ceiling: "
"this is SCALAR (Nordstrom) gravity — it gives Newton's law but NOT light bending, and not two-polarization "
"gravitational waves. Full general relativity requires a massless SPIN-2 field protected by DIFFEOMORPHISM "
"INVARIANCE (Weinberg's uniqueness theorem), and a fixed-background medium has no diffeomorphism invariance — so GR "
"itself remains out of reach. It also requires the medium to sit near criticality, a fine-tuning — though that is "
"arguably the emergent statement of WHY GRAVITY IS SO WEAK, and it is a prediction rather than a fudge.")

heading("8.11  From scalar to tensor gravity: two gravities, and deconfinement", 2)
body("Section 8.10 gives a SCALAR gravity, and the honest next question is whether the medium can reach the TENSOR "
"theory — general relativity, with light bending and spin-2 waves. The model in fact carries TWO universal "
"attractions that disagree on the light-bending parameter gamma: the amplitude mode of 8.10 (spin-0, gamma = 0) and "
"the medium's incompatible-strain / curvature sector (spin-2, gamma = 1). Which one dominates the long-range force is "
"decided not by the couplings but by the MASS: a massless mediator's 1/r always beats a massive mediator's Yukawa. So "
"if the graviton is massless and the amplitude mode is gapped, gamma climbs from 0 to 1 across the amplitude Compton "
"wavelength 1/m_A — a falsifiable SCALE-DEPENDENT gamma, GR at long range and scalar-contaminated below 1/m_A. This "
"also reinterprets 8.10: making the scalar long-range by tuning m_A -> 0 forces gamma = 1/2, which is observationally "
"excluded — so the amplitude mode should stay gapped, and gravity proper is the spin-2 graviton.")
body("But the spin-2 sector is CONFINING in the pure medium. Its energy is the elastic cost of incompatible strain — "
"the biharmonic (Kirchhoff-plate) action kappa (lap phi)^2 — whose 3D Green's function grows LINEARLY with "
"separation, giving a constant force between two curvature charges: a string tension, not a Newtonian force. To become "
"gravity it must DECONFINE. The mechanism is the same Sakharov loop that induced the photon's Maxwell term (8.7): "
"integrating out the gapped matter generates an Einstein-Hilbert term, which at quadratic order adds an ordinary "
"two-derivative stiffness mu (grad phi)^2 to the biharmonic. The full propagator is then 1/(kappa q^4 + mu q^2), "
"whose exact three-dimensional Green's function is closed-form:")
add_eq("G(R) = (1 / 4 pi mu R)(1 - exp(-R / ell)),     ell = sqrt(kappa / mu)", "8.11a")
body("The lower-derivative induced term dominates the infrared: below the crossover ell the force is the confining "
"string tension, and ABOVE ell it turns over into an exact inverse-square Newtonian tail with G = 1/(4 pi mu). "
"Measured on the correct tool for a spherically symmetric source — a one-dimensional radial ODE with no box and no "
"periodic images — the Newton-tail exponent is -2.0000, Newton's constant matches 1/(4 pi mu) to five figures, the "
"closed form (8.11a) is reproduced to one part in ten million, and the pure-medium (mu = 0) confinement shows its "
"linear growth directly, G ~ R^+1.000 — the clean form of the string tension that no finite box can display. ANY "
"positive mu deconfines the graviton: the confining +R becomes a Newtonian -1/r.")
result("Result 8.11 — the graviton deconfines.", "The tensor (curvature) sector is confining in the pure medium — a "
"constant-force string tension, the clean 3D form of the earlier 'curvature charges repel and grow' result. A "
"positive Sakharov-induced Einstein term turns it into a massless Newtonian graviton: the force crosses over at "
"ell = sqrt(kappa/mu) from the string tension to an exact 1/r^2 with G = 1/(4 pi mu). Everything now rests on the SIGN "
"of the induced mu.")

heading("8.12  The sign of induced gravity, and a dynamical spin-2 graviton", 2)
body("Deconfinement needs mu > 0, and this is the notorious Sakharov sign — the induced Einstein coefficient is "
"ultraviolet-dominated and scheme-sensitive, and free fields do not universally give the healthy sign. The model "
"settles it by a calibration it already owns: its induced PHOTON is healthy (8.7), so the induced Coulomb kinetic "
"term — the charge-density correlator <J0 J0> — is a healthy dielectric, susceptibility chi > 0. Computing the "
"induced NEWTONIAN kinetic term — the energy-density correlator <T00 T00>, since energy is the gravitational charge — "
"from the SAME gapped-Dirac loop with identical conventions, its momentum-squared coefficient comes out the SAME "
"(positive) sign as the photon's in every case tested: five masses, three cutoffs, and both energy-density vertex "
"definitions. Induced gravity is therefore as healthy as the electromagnetism the model already runs on: mu > 0.")
body("What kind of graviton is it? On the 2+1D emergent cone the SPATIAL spin-2 graviton is non-dynamical — a massless "
"symmetric tensor has D(D-3)/2 physical polarizations, which is zero in three spacetime dimensions and two in four. "
"So the radiative, light-bending part of gravity can only appear in 3+1D, and there it does: a four-component Dirac "
"loop induces a NONZERO transverse-traceless kinetic term (whereas the same construction gives zero in 2+1D), its two "
"polarizations h_plus and h_cross are degenerate to four figures (a single helicity-2 field, not two unrelated "
"modes), and its sign matches the induced transverse photon (a healthy, not ghost, mode). Both sectors of the "
"graviton — the Newtonian h00 and the radiative spin-2 — are thus induced and healthy in the physical dimension.")
body("The Einstein normalization gamma = 1 is the transversality (Ward identity) of the induced graviton, and here the "
"model draws a sharp and honest line. On a periodic lattice (a torus, with no boundary and hence no cutoff surface "
"term) the U(1) photon Ward identity closes to MACHINE PRECISION once the diamagnetic seagull is included — the "
"symmetry-preserving regulator works, because U(1) gauge invariance is an EXACT lattice symmetry. Diffeomorphism "
"invariance is NOT a lattice symmetry — a lattice keeps only discrete translations — so the graviton Ward identity "
"cannot close exactly: it is inhomogeneous (the induced vacuum stress <T> is nonzero, the same term that reappears as "
"the cosmological constant in 8.13), and no finite-cutoff DIRECT measurement can force gamma = 1. That is not a "
"numerical shortcoming but a structural fact: gamma = 1 is EMERGENT in the infrared, on the same footing as the "
"model's emergent Lorentz invariance. It is read from Weinberg's theorem — a massless spin-2 coupled to the CONSERVED "
"infrared (Dirac) stress tensor is forced to be Einstein — now with the previously-missing ingredient supplied by "
"measurement: the spin-2 graviton genuinely propagates and is healthy.")
result("Result 8.12 — induced gravity is healthy, spin-2, and Einstein in the infrared.", "The induced Newtonian "
"coupling has the same sign as the model's working photon: mu > 0, so the graviton of 8.11 deconfines into real "
"attraction. In 3+1D the radiative spin-2 graviton is dynamical, doubly degenerate (helicity 2), and healthy. And "
"gamma = 1 follows from Weinberg on the conserved infrared stress tensor — an EMERGENT identity, since diffeomorphism "
"invariance (unlike the exactly-closing U(1) Ward identity) is not a lattice symmetry. General relativity is thus "
"reached as an infrared fixed point. Open still: the MAGNITUDE of G is cutoff-dependent (the Sakharov feature), and a "
"lattice-exact gamma = 1 is structurally unavailable, by design.")

heading("8.13  The cosmological constant", 2)
body("The same induced vacuum stress that made the graviton Ward identity inhomogeneous is, physically, the "
"cosmological constant — and it is the sharpest quantitative disaster in physics. The medium's zero-point energy is "
"of order the microscopic (node) scale; with the node spacing fixed at the Planck length (8.8's scale-fixing), that "
"is about 10^122 times the observed dark-energy density. Taken at face value it demands a 122-digit fine-tuning. But "
"the model's vacuum is not empty space with fields on top — it is a self-sustained CONDENSATE, and that changes what "
"gravitates. Following Volovik's analysis of emergent gravity in quantum liquids, the emergent metric couples to the "
"vacuum STRESS — the grand-canonical potential density rho_Lambda = eps - mu n = -P — not to the bare energy density "
"eps. A self-sustained vacuum, one that can exist with nothing outside pushing on it (which is what the vacuum of "
"empty space is), has zero pressure:")
add_eq("rho_Lambda = eps - mu n = -P;      self-sustained vacuum: P = 0  =>  rho_Lambda = 0", "8.13a")
body("The huge zero-point energy is absorbed into the equilibrium condensate density, not into curvature. This is not "
"a tuning: the density SELF-ADJUSTS so that P = 0, for ANY bare eps. Sweeping the bare vacuum energy across all 122 "
"orders of magnitude, the gravitating rho_Lambda stays zero to machine precision while a rigid (non-adjusting) vacuum "
"would gravitate the full eps — the standard disaster. The observed small but nonzero Lambda is then the residual of "
"a slight DEPARTURE from equilibrium (rho_Lambda scales with that departure), not a cancellation of 122 digits.")
result("Result 8.13 — the cosmological-constant fine-tuning, dissolved (not the value derived).", "Because the vacuum "
"is a self-sustained condensate, the quantity that gravitates is its grand potential -P, which vanishes at "
"equilibrium — automatically, for any bare zero-point energy, with the density self-adjusting. The 10^122 fine-tuning "
"is thus dissolved: the equilibrium vacuum gravitates NOTHING by thermodynamics. Honest ceiling: this does not "
"predict the observed nonzero Lambda, which is relocated to a cosmological question — why the vacuum sits slightly "
"off equilibrium (expansion, matter, relaxation) — and remains open.")

heading("8.14  The Standard-Model gauge group: emergent Yang-Mills", 2)
body("The emergent photon of Section 8.7 is Abelian. The Standard Model needs the NON-ABELIAN groups SU(2) and SU(3), "
"and the question is whether the same fermion loop gives genuine Yang-Mills or merely several decoupled copies of the "
"photon. The decisive difference is the gauge-boson SELF-INTERACTION, and it has a clean signature. A non-Abelian "
"field strength carries a commutator, F = dA + i[A, A], so a spatially UNIFORM non-Abelian field has nonzero field "
"strength from that commutator alone, while a uniform Abelian field is always pure gauge. Placing a fermion in the "
"fundamental of SU(N) with uniform background links and measuring the induced action (the filled-sea energy), a "
"COMMUTING (Cartan) configuration costs nothing — pure gauge, to machine precision at every amplitude — while a "
"NON-COMMUTING configuration costs an induced action that grows as A^4, exactly the Yang-Mills Tr[A_x, A_y]^2. N^2-1 "
"decoupled photons would give zero for BOTH; the A^4 cost is the fingerprint of a genuine, self-interacting "
"non-Abelian field. It appears for SU(2) (three gauge bosons) and SU(3) (eight gluons), with a single universal "
"coupling guaranteed by the exact non-Abelian lattice gauge invariance of the Wilson links — the same exact-symmetry "
"footing that made the U(1) Ward identity close exactly in 8.12.")
result("Result 8.14 — emergent non-Abelian gauge fields.", "The Sakharov mechanism that induced the photon induces "
"genuine Yang-Mills: a uniform non-commuting gauge field costs an induced action ~ A^4 = Tr[A,A]^2 (the "
"self-interaction), while a commuting one is exactly pure gauge — shown for SU(2) and SU(3), with a universal coupling "
"from exact lattice gauge invariance. Emergent gauge theory thus scales from U(1) to SU(N). Honest ceiling: this is "
"the MECHANISM, not the Standard Model — it does not derive the specific group SU(3)xSU(2)xU(1), the chiral coupling, "
"the fermion representations or hypercharges, anomaly cancellation, or the Higgs. The group is an input; its "
"Yang-Mills dynamics, and the fermion light cone they inherit, are induced.")

heading("8.15  Chirality without inconsistency: the anomaly and its inflow", 2)
body("Section 8.14 supplies non-Abelian gauge fields and Section 8.2 a single chiral fermion on a domain wall. The "
"Standard Model needs BOTH AT ONCE, and that is where the real obstruction sits. The SM's SU(2) couples only to "
"LEFT-handed fermions, and a chiral gauge theory is not merely incomplete but INCONSISTENT unless its anomalies "
"cancel: a chiral fermion's gauge current is not conserved in a background field, so a lone chiral fermion coupled to "
"a gauge field is not a theory at all. The domain-wall construction survives this by a specific mechanism -- "
"CALLAN-HARVEY ANOMALY INFLOW. Each wall is INDIVIDUALLY anomalous, and the charge it appears to lose is supplied by "
"the BULK, which pumps it to the other wall. The content is quantitative and quantized: the bulk Chern number, the "
"number of chiral modes per wall, and the charge pumped per flux quantum are ONE integer.")
body("Measured on the same Wilson-Dirac strip that carried the chiral fermion: in the topological phase the bulk Chern "
"number is -1.000, the edge spectral flow is +1 on one wall and -1 on the other, and the two SUM TO ZERO; in the "
"trivial phase all three vanish. So each wall separately violates charge conservation by exactly one unit per flux "
"quantum -- a genuine anomaly -- while the lattice as a whole is vector-like and anomaly-free, exactly as "
"Nielsen-Ninomiya demands. Neither wall is a consistent theory alone; the wall pair together with the bulk is.")
result("Result 8.15 — chirality is consistently realizable, by inflow.", "The anomaly here is not a pathology but a "
"quantized bookkeeping identity, and the three numbers that must agree do agree. This is the mechanism "
"Standard-Model chirality requires, and it operates in the model's own domain-wall construction alongside the induced "
"non-Abelian gauge fields of 8.14. HONEST ceiling: this is NOT the Standard Model's anomaly cancellation. The SM is a "
"standalone FOUR-dimensional chiral gauge theory whose anomalies cancel among its OWN fermion content -- the "
"quark/lepton hypercharge conspiracy -- with no bulk to lean on. Here the BULK does the cancelling, so the wall "
"theory is anomaly-free only together with it. Producing a standalone anomaly-free chiral spectrum, i.e. the actual "
"SM fermion content and hypercharges, is not attempted and is not fixed by the medium.")

heading("8.16  The model against the data, and a retraction", 2)
body("A model earns its keep by being falsifiable, and the two predictions this project carries had never been put "
"against real numbers. Doing so sharpens one, deflates the other, and forces a retraction.")
body("The LORENTZ-VIOLATION signature of Section 8.8 is QUADRATIC and Planck-suppressed. Compared with photon "
"time-of-flight limits on quadratic dispersion, which the literature places at an effective scale of order 1e10 to "
"1e12 GeV, the model's effect is about 1e-16 of current sensitivity. It is therefore safe against every existing "
"bound by roughly sixteen orders of magnitude -- but that cuts both ways, and the honest conclusion is that this "
"prediction is NOT currently falsifiable: reaching it would need some eight orders of improvement in quadratic-"
"dispersion sensitivity. Present experiments are sensitive to a LINEAR, n = 1 signal, which this model specifically "
"does not predict. Section 8.8's framing of this as the project's testable prediction was too generous.")
body("The GRAVITATIONAL prediction is the real one. With a massless graviton and a gapped amplitude mode (8.11), the "
"amplitude mode survives as a Yukawa ADDITION to gravity of range 1/m_A, so gamma falls below 1 inside that range. "
"Cassini's solar-system bound turns out to be nearly useless here, permitting lambda up to about 1e10 m; the binding "
"constraint is short-range gravity, where torsion-balance tests of the inverse-square law beat Cassini by fourteen "
"orders of magnitude. Taking the Yukawa strength of order unity, as is natural when both mediators couple to energy, "
"the amplitude gap must satisfy m_A of order 4 meV or more -- which sits within a factor of about 1.6 of the "
"dark-energy scale, 2.4 meV, the well-known coincidence that makes the sub-millimetre range the frontier of these "
"experiments. This is a genuine, reachable test: a gravitational-strength Yukawa just below current reach would show "
"up as gamma < 1 at short distance.")
body("The same confrontation RETRACTS the 1e122 criticality tuning reported earlier in this program. That number "
"followed from reading gravity's RANGE as 1/m_A -- the scalar picture -- so that any bound on gravity remaining "
"inverse-square out to large distances forced the amplitude gap to be absurdly small. The tensor arc of Sections "
"8.11 and 8.12 replaced that premise: the long-range force is carried by the MASSLESS deconfined graviton, and the "
"amplitude mode is only a short-range correction. The surviving experimental constraint is therefore a LOWER bound "
"on m_A, with no upper bound at all -- a larger gap is only safer -- and an UNTUNED medium clears it by some thirty "
"orders of magnitude. The fine-tuning is dissolved, not reduced. The scale-fixing result that a0 = l_Planck is "
"unaffected and stands.")
result("Result 8.16 — one prediction sharpened, one deflated, one tuning retracted.", "Confronted with data: the "
"quadratic Lorentz violation is safe by ~16 orders and, honestly, out of experimental reach -- not the testable "
"prediction it was billed as. The gravitational prediction IS testable: short-range tests already require the "
"amplitude gap to exceed about 4 meV, coincidentally within a factor 1.6 of the dark-energy scale and squarely in "
"the sub-millimetre window now being probed. And the 1e122 criticality tuning is RETRACTED, because it rested on the "
"superseded scalar reading of gravity. Together with Section 8.13, which dissolved the cosmological-constant tuning "
"by the condensate's equilibrium thermodynamics, the model has now shed BOTH of its 10^122 fine-tunings -- each "
"removed by a structural result rather than a fitted parameter. (Scope: literature order-of-magnitude bounds, "
"assuming comparable scalar and tensor couplings; loosening that weakens the gravitational bound proportionally.)")

heading("8.17  The first dynamical integration: chiral matter and a gauge field, in real time", 2)
body("Every result to this point -- emergent Lorentz invariance, chiral fermions, the induced photon and graviton, the "
"anomaly -- is established at the level of a dispersion relation, a band structure, or a defect algebra. None is a "
"RUNNING SIMULATION in which two of these emergent sectors coexist and interact IN TIME. That gap is the sharpest "
"honest criticism of the whole program, and this section takes the first step across it, using the anomaly of "
"Section 8.15 as the target. That anomaly was established statically, by counting; the question here is whether the "
"charge it bookkeeps actually MOVES.")
body("The test is a Laughlin flux threading. On the same QWZ strip -- two chiral walls, x periodic -- the "
"negative-energy sea is filled at t = 0, and one full flux quantum is threaded adiabatically by ramping a uniform "
"vector potential, A: 0 -> 2 pi / L_x, which enters as k_x -> k_x + A(t). Every occupied orbital is evolved under the "
"time-dependent Schrodinger equation by exact exponentiation of the instantaneous Hamiltonian at each step -- no "
"adiabatic-following shortcut, which would presuppose the result. If the anomaly is genuine dynamics rather than "
"bookkeeping, exactly one unit of charge must cross from one wall to the other per flux quantum, pumped through the "
"BULK since the walls are spatially separated and nothing local joins them.")
body("It does. The charge in the bottom half of the strip changes by 0.9991 of one unit in the topological phase and "
"by exactly zero in the trivial control, so the transfer is the anomaly and not the ramp. Slowing the ramp drives the "
"residual to zero as roughly one over the number of steps -- the quantization is physics, not an artifact of the "
"discretisation. The charge has crossed the bulk in real time: Callan-Harvey inflow, observed as dynamics.")
result("Result 8.17 — the anomaly happens, and two sectors run together.", "Threading one flux quantum through the "
"chiral strip and evolving the filled sea under the actual time-dependent Schrodinger equation pumps exactly one unit "
"of charge between the walls (0.9991, converging to one as the ramp slows), and exactly zero in the trivial control. "
"The topological accounting of Section 8.15 therefore describes real time evolution -- the bulk really does supply "
"what each wall loses. This is the program's FIRST running simulation in which two emergent sectors -- chiral matter "
"and a gauge field -- coexist and interact in time, and the consistency does not break down. HONEST scope: the field "
"threaded is the U(1) GAUGE field, not gravity. The full integration this program's limitations call for -- emergent "
"Lorentz-invariant chiral quantum matter interacting through an emergent SPIN-2 gravity, with back-reaction -- is not "
"done; gravitational back-reaction in a running simulation remains the open integration problem. What is shown is that "
"the first pair of sectors can be run together at all.")

heading("8.18  Gravitational back-reaction, run as a conserving simulation", 2)
body("Section 8.17 ran chiral matter together with a gauge field. GRAVITY is the harder and more important case, "
"because back-reaction IS gravity's defining feature: matter tells geometry how to curve, geometry tells matter how to "
"move, and the two must be solved TOGETHER and self-consistently. Until that is done in time, a gravity result is a "
"dispersion relation rather than a force. This section does it, and it is the step at which the gravity sector stops "
"being a toy.")
body("The system evolved is the Schrodinger-Newton pair -- the non-relativistic limit of a massive matter field "
"minimally coupled to its own gravity, and the standard model of self-gravitating quantum matter:")
add_eq("i d_{t} psi = -(1/2) lap psi + Phi psi,      lap Phi = 4 pi G |psi|^{2}", "8.18a")
body("The gravity in (8.18a) is not inserted by hand: it is the infrared-effective form of the gravity this program "
"DERIVED -- the deconfined graviton mediates an exact Newtonian potential with G = 1/(4 pi mu) (Section 8.11), and the "
"sign mu > 0 was measured against the model's own healthy photon (Section 8.12). Co-evolving matter with that "
"potential uses the derived gravity rather than inventing one.")
body("Four properties separate a scientific simulation from a toy, and all four hold. CONSERVATION: evolved by a "
"split-step scheme that is symplectic in the matter sector with the potential resolved self-consistently each step, "
"the total energy is conserved to three parts in a billion and the norm to four parts in a hundred trillion, through "
"the full nonlinear evolution -- against which the earlier gravity-by-density attempt, which leaked some eighty per "
"cent of its energy through an inconsistent cutoff, is the instructive failure. SELF-BINDING: with gravity on, a "
"packet forms a self-gravitating BOUND STATE of negative total energy whose width settles about a finite soliton "
"scale, while the identical packet with gravity off disperses without bound -- the single self-bound lump that the "
"nearly incompressible medium of Phase 3 could never produce. EQUILIBRIUM: imaginary-time relaxation finds the soliton "
"ground state, and it satisfies the scale-virial identity 2T + W = 0 to within three hundredths, the signature of a "
"genuine gravitational equilibrium rather than a long-lived transient. This last check requires ISOLATED, free-space "
"boundary conditions; a periodic box distorts the long-range potential and leaves a spurious virial residual near one, "
"which is reported alongside as a methodological control. CONVERGENCE: the relaxed soliton's energy and virial settle "
"monotonically as the mesh is refined, so the bound state is a property of the continuum system and not of the grid.")
result("Result 8.18 — gravity as a force that can be run.", "The coupled matter-plus-gravity system evolves "
"self-consistently with energy conserved to ~1e-9 and norm to ~1e-14; it BINDS matter into a self-gravitating soliton "
"that a gravity-off control does not form; the relaxed soliton satisfies the virial identity 2T + W = 0; and both "
"converge under mesh refinement. This is the gravity the program derived, run as a FORCE in time rather than read off "
"a propagator -- the first back-reacting gravitational simulation here, and the point at which the gravity sector "
"meets the standards of a scientific simulation. HONEST scope: non-relativistic (Schrodinger, not Dirac matter), "
"scalar/Newtonian (the h00 sector; the radiative spin-2 graviton of Section 8.12 is not evolved), and the matter is a "
"classical field. Chiral QUANTUM matter interacting through the emergent SPIN-2 gravity with radiative back-reaction "
"is still the open integration problem.")

heading("8.19  The magnitude of Newton's constant", 2)
body("One caveat has trailed every gravity result: the STRENGTH of gravity was called cutoff-dependent. That is the "
"standard Sakharov ambiguity -- the induced Einstein-Hilbert coefficient is ultraviolet-dominated, so in a continuum "
"theory with an arbitrary cutoff its magnitude is arbitrary. The caveat does not apply to this model, because the "
"medium has a PHYSICAL ultraviolet cutoff: the node spacing a0, fixed to the Planck length in Section 8.8. Newton's "
"constant is therefore not a free parameter but a definite number, computable from the lattice.")
body("Gravity's mediator takes its kinetic term from the fermion loop, and G is set by the induced Newtonian stiffness "
"mu -- the coefficient of q^2 in the energy-density correlator, the h00 sector's induced 1/(4 pi G). Evaluated over "
"the WHOLE Brillouin zone, so that the cutoff is the physical lattice scale rather than an arbitrary disc, the loop "
"returns a number of order unity in lattice units. Hence G = 1/(4 pi mu) is of order a0^2, that is of order the Planck "
"area, up to an order-unity factor -- computed rather than fitted. Including more light fermion species stiffens the "
"geometry in the standard way, G ~ a0^2 / N_f, so even a Standard-Model-like species count leaves G at a few per cent "
"of the Planck area: still Planckian. Gravity is weak for exactly one reason, that a0 is Planckian, and the induced "
"coefficient supplies no hierarchy and requires no tuning.")
result("Result 8.19 — G is the Planck area, with no tuning.", "With a physical cutoff the induced stiffness is an "
"order-unity number and G = O(1) x a0^2 = O(1) x l_Planck^2. The weakness of gravity is entirely the smallness of the "
"Planck length; the species count only sharpens it. Taken with Section 8.8, which fixed a0 = l_Planck by matching the "
"measured G, the loop closes: the node spacing IS the Planck length, and the gravity the medium induces has Planck "
"strength. HONEST scope: the SCALE (order a0^2) and the SIGN (positive, hence attractive) are robust, but the precise "
"order-unity coefficient remains scheme-sensitive -- the residual of the Sakharov ambiguity, since the lattice is one "
"regulator among several. What is settled is the qualitative point that had been left open: the magnitude of G is not "
"free here.")

heading("8.20  Gravity that radiates, and the monopole that cannot", 2)
body("Section 8.18 ran gravity as a force, but in the NEWTONIAN limit: the potential was solved from an instantaneous "
"constraint, which is exact only for slow sources. A real gravitational field is RETARDED -- it propagates at finite "
"speed and carries energy away as radiation -- and that radiative sector is the last structural piece of the gravity "
"programme. It also carries a sharp signature which separates the TENSOR gravity this programme arrived at from the "
"SCALAR gravity it discarded along the way.")
body("The signature is the MONOPOLE. A scalar (Nordstrom) gravity radiates from a spherically pulsating mass: a "
"breathing star emits scalar gravitational waves. A spin-2 field cannot. The radiative degrees of freedom of "
"linearised gravity are the transverse-traceless part of h_ij, and the TT projection annihilates a spherically "
"symmetric source identically, so monopole radiation is forbidden and the leading channel is the quadrupole. That is "
"why a pulsating star does not gravitationally radiate, and it is a structural test of whether the model's gravity is "
"really spin-2.")
body("Evolving the linearised TT wave equation in momentum space -- each mode a driven oscillator, the source switched "
"off after a few cycles, and the field energy remaining once the near field has dispersed counted as the radiated "
"energy -- gives three results. The disturbance PROPAGATES: after the source stops, the outgoing shell moves at speed "
"0.96 in units where c = 1, the few-per-cent deficit being the finite width of the shell rather than a slow wave. The "
"radiation is carried by exactly TWO polarizations, the transverse-traceless projector having rank two, which are the "
"helicity-2 states of Section 8.12 and not a scalar breathing mode. And at identical amplitude, width and frequency, "
"the monopole source radiates about 1e-14 of what the quadrupole source radiates -- a ratio of 3e-13, which is machine "
"zero against a finite signal.")
result("Result 8.20 — the radiation is spin-2, and the monopole channel is closed.", "The model's gravitational field "
"propagates at c and carries exactly two polarizations, and a spherically pulsating mass radiates NOTHING through it, "
"while an equal quadrupole radiates a finite amount. Scalar gravity would have opened the monopole channel; spin-2 "
"gravity forbids it, and the measurement finds it shut to machine precision. This is a structural confirmation, "
"independent of the earlier propagator arguments, that the gravity this programme ended up with is genuinely tensor. "
"HONEST scope: the radiation is LINEARISED and the source is PRESCRIBED, so the back-reaction of the radiated energy "
"ON the source -- the inspiral of a binary -- is not computed. Neither is the quadrupole LUMINOSITY formula tested: "
"that law assumes a source small compared with the wavelength, the source used here is not, and the frequency "
"dependence in this setup therefore reflects the source's own spatial spectrum rather than the multipole expansion. No "
"such claim is made.")

result("Result 8 — scorecard.", "The barriers usually fatal to a 'space is a medium' theory now carry concrete "
"in-model demonstrations: emergent Lorentz invariance, emergent fermions (a Dirac cone plus a single chiral "
"fermion on a domain wall), a proper relativistic QFT on quantization, and an emergent photon. More striking than "
"the individual results is that ONE principle — everything from one structure — surmounts several at once. Even quantum "
"mechanics is largely condensate mechanics (§8.6): the Schrodinger wave, hbar, the Born rule and de Broglie's "
"lambda = h/p emerge, leaving only the measurement problem's hard core as a postulate. Cast against experiment, the "
"model yields a specific, falsifiable Lorentz-violation signature (§8.8), its first genuine prediction. GRAVITY "
"took the longest and cost a retraction: §8.9 records the failure of every ELASTIC route (the curvature sector is "
"unshieldable but REPULSIVE and growing; the tetrad is long-range in FIELD but shieldable and force-free by "
"Eshelby-Crum), and §8.10 then SOLVES it by applying the project's own range principle to gravity's mediator for the "
"first time — the condensate's unprotected, hence gapped, AMPLITUDE mode, which is exactly why gravity always looked "
"screened, and which at criticality yields Newton's law, 1/r^2, universally attractive. Sections 8.11-8.12 then carry "
"gravity from that SCALAR force to the TENSOR theory: the confining curvature sector DECONFINES into a Newtonian "
"graviton once the Sakharov loop supplies a positive Einstein term (and that sign is MEASURED, mu > 0, by calibration "
"against the model's own healthy photon), the radiative spin-2 graviton is dynamical and healthy in 3+1D, and gamma = 1 "
"follows from Weinberg as an INFRARED-emergent identity — general relativity as a fixed point rather than a "
"lattice-exact law. Section 8.13 dissolves the cosmological-constant fine-tuning (the self-sustained condensate "
"vacuum gravitates its grand potential -P, which vanishes at equilibrium for any zero-point energy), and Section 8.14 "
"shows the photon's induction mechanism scales to non-Abelian YANG-MILLS (SU(2), SU(3)), while Section 8.15 shows "
"CHIRALITY is consistently realizable alongside it, the anomaly being a quantized identity settled by Callan-Harvey "
"inflow -- now confirmed DYNAMICALLY (Section 8.17), in the program's first running simulation of two emergent sectors "
"together. Section 8.16 then puts the model against data, and the result is bracing: the Lorentz-violation signature is "
"safe but NOT presently falsifiable, the genuinely testable prediction is the short-range gravitational gamma, and "
"the 1e122 criticality tuning is RETRACTED. Sections 8.18-8.19 then close two more: gravitational BACK-REACTION now runs as a conserving, convergent simulation that binds matter into a virial-satisfying soliton, and the magnitude of G is fixed at O(1) x the Planck area by the physical lattice cutoff. Section 8.20 then shows the field RADIATES like spin-2 -- propagating at c with two polarizations, and with the monopole channel shut to machine precision, which scalar gravity would have left open. What remains open is the DERIVATION of the specific Standard-Model group, representations, hypercharges and chiral "
"content (all still inputs), the observed value of the cosmological constant, and the measurement problem's hard "
"core. Four self-corrections are now on record — the within-sector Lorentz result, the retracted gravity route, the "
"refuted critical-nucleus prediction, and the retracted 1e122 tuning — which is the discipline working as intended.")
table(["Barrier","Status","Key result"],
 [["Emergent Lorentz","achieved","one round universal cone, violations ~ (E/E_Planck)^2; cross-statistics universality holds once the boson is a fermion composite (§8.5)"],
  ["Fermions","achieved","Dirac cone on honeycomb; single chiral fermion on a domain wall (evades Nielsen-Ninomiya)"],
  ["Quantum mechanics","largely emergent","quantizes to a relativistic QFT; and from the condensate directly (§8.6): the Schrodinger wave + hbar as a material property, the Born rule as a stochastic attractor, and de Broglie v=grad(S)/m for a particle's own wave. Only guidance by a separate pilot wave + definite outcomes (measurement) stay a postulate"],
  ["Emergent photon","achieved","the Dirac-node position: a fluctuation of the medium's own bonds, on the fermion cone"],
  ["Long-range gravity","achieved as SCALAR gravity (§8.10)","the mediator is the condensate's AMPLITUDE mode: unprotected, hence gapped — which is WHY gravity always looked screened. It couples monopolarly to positive-definite energy, and scalar exchange between like charges ATTRACTS. Measured: lambda*m_A = 1.00, and E = -C exp(-R/lambda)/R exactly, so at criticality Newton's law, 1/r^2, universally attractive. The ELASTIC route is provably dead (§8.9: Bitter-Crum + Eshelby-Crum)"],
  ["Full general relativity","reached as an IR fixed point (§8.11-8.12)","the confining curvature sector DECONFINES into a Newtonian graviton once the induced Einstein term mu>0 (measured, by calibration against the healthy photon); gamma=1 follows from Weinberg on the conserved IR stress tensor. Diffeomorphism invariance is EMERGENT (not a lattice symmetry), so gamma=1 is an IR identity, not lattice-exact; the magnitude of G stays cutoff-dependent"],
  ["Spin-2 graviton (dynamical)","achieved in 3+1D (§8.12)","the transverse-traceless graviton is NON-dynamical in 2+1D (0 polarizations) but dynamical in 3+1D (2 polarizations): the induced TT kinetic term is nonzero, the two polarizations are degenerate (helicity 2), and the mode is healthy (same sign as the transverse photon)"],
  ["Empirical prediction","one testable, one out of reach (§8.16)","the n=2 Lorentz violation is safe vs every bound by ~16 orders — but being QUADRATIC it is ~1e-16 below current sensitivity and NOT presently falsifiable (§8.8 tempered). The reachable prediction is gravitational: a scale-dependent gamma below 1/m_A (§8.11), for which short-range tests already require m_A >~ 4 meV — within a factor 1.6 of the dark-energy scale, in the sub-mm window now being probed"],
  ["Gravitational back-reaction","runs as a conserving simulation (§8.18)","matter sources the potential and the potential moves matter, solved together: energy conserved to ~1e-9 and norm to ~1e-14, a self-gravitating BOUND soliton forms (a gravity-off control disperses), the relaxed soliton satisfies the virial identity 2T+W=0, and both converge under mesh refinement. Scope: non-relativistic, scalar/Newtonian, classical matter; the radiative spin-2 sector is not evolved"],
  ["Magnitude of G","fixed at the Planck area (§8.19)","the Sakharov cutoff-ambiguity does not apply because the cutoff is PHYSICAL (a0 = l_Planck): over the full Brillouin zone the induced stiffness is O(1) in lattice units, so G = O(1) a0^2, with G ~ a0^2/N_f. Gravity is weak because a0 is Planckian -- no hierarchy, no tuning. The O(1) coefficient stays scheme-sensitive"],
  ["Gravitational radiation","spin-2, monopole forbidden (§8.20)","the linearised TT field propagates at c with exactly 2 polarizations, and a spherically pulsating source radiates ~1e-14 of an equal quadrupole (ratio 3e-13, machine zero). Scalar gravity would radiate the monopole; spin-2 forbids it. Scope: linearised, prescribed source -- inspiral back-reaction and the quadrupole luminosity law are NOT tested"],
  ["Chirality + anomalies","consistent by inflow (§8.15)","a chiral gauge theory is inconsistent unless anomalies cancel; here bulk Chern number = chiral modes per wall = charge pumped per flux quantum = one integer (measured -1, +1/-1, sum 0). Each wall is anomalous, the lattice is vector-like, the bulk supplies the inflow. NOT the SM's own 4D cancellation — the bulk does the work"],
  ["Cosmological constant","fine-tuning dissolved (§8.13)","the self-sustained condensate vacuum gravitates its grand potential -P, which vanishes at equilibrium for ANY bare zero-point energy (measured across 122 orders, no tuning). The equilibrium value is exactly zero; the observed nonzero Lambda is relocated to a departure from equilibrium (open)"],
  ["Non-Abelian gauge fields","mechanism achieved (§8.14)","the fermion loop induces genuine Yang-Mills for SU(2) and SU(3): a uniform non-commuting field costs ~A^4 = Tr[A,A]^2 (self-interaction), a commuting one is pure gauge. Universal coupling from exact lattice gauge invariance"],
  ["SM group / reps / constants","open","the GROUP is an input; the derivation of SU(3)xSU(2)xU(1), the chiral coupling, fermion representations, hypercharges, anomaly cancellation, and the Higgs is not attempted"]],
 cap="Table 5.  The fundamental-physics barriers and their status in the model.")

# ===== 9 Synthesis =====
heading("9  Synthesis", 1)
body("Read as a whole, the program builds one object — a medium carrying a field — and promotes it step by "
"step: a real field (particles), a complex field (charge and spin), a mobile medium (robust topology), a "
"density coupling (gravity), a larger target space (3D point charges), and a gauge field (electromagnetism). "
"At each step the same question recurs and receives the same kind of answer. A bare topological charge's "
"self-energy always grows with system size — logarithmically for a 2D vortex, linearly for a 3D line or a 3D "
"hedgehog — and a strong long-range force between charges appears only when the mediator is massless and "
"protected by a symmetry: a global Goldstone mode, or a gauge photon in the Coulomb phase. Gauge the same "
"symmetry into its broken phase and the force screens. From one medium and one field, the model thus "
"reproduces the qualitative menu of real forces and both gauge phases, with the reach of each governed by a "
"single principle.")

# ===== 10 Limitations =====
heading("10  Limitations and scope", 1)
body("This is a toy model, and the results are qualitative correspondences, not derivations of the Standard "
"Model or general relativity. The gravity sector is weak and non-self-binding on the nearly incompressible "
"medium; a definitive compressible-medium test awaits a properly conservative soft potential. The 3D "
"point-charge interaction of Route B is measured on a pinned, relaxed texture, and Route C's Coulomb phase is "
"the classical minimum of the Maxwell energy rather than a full quantum treatment. Range claims of the "
"long-range kind rest on system-size scaling rather than infinite volume. The correspondences are "
"qualitative — 'EM-like', 'nuclear-like' — and quantitative matching of coupling constants or mass ratios is "
"not attempted. What survives these caveats is the central, dimension-robust and repeatedly gate-checked "
"statement: within this medium, a force's range is set by whether a symmetry protects its mediator from a "
"mass term.")
body("The fundamental-physics program of Section 8 carries a further, honest qualification. Most of those results are "
"established at the level of the dispersion relation, the band structure, and the defect algebra rather than a single "
"running simulation. Sections 8.17 and 8.18 take the first steps past this. Section 8.17 runs chiral matter and a "
"gauge field together in time and recovers the anomaly as real charge transport; Section 8.18 runs GRAVITY with "
"genuine back-reaction -- matter sourcing the potential and the potential moving matter, self-consistently, with "
"energy conserved to parts in a billion, a self-gravitating bound state formed, the virial identity satisfied, and "
"convergence under mesh refinement. Those are the standards a scientific simulation must meet, and in those sectors "
"the model now meets them. What is still NOT a single running simulation is the full target: emergent "
"Lorentz-invariant, chiral, QUANTUM matter interacting through the emergent SPIN-2 gravity with radiative "
"back-reaction. Section 8.18 is non-relativistic, scalar/Newtonian and classical in its matter; Section 8.20 evolves "
"the radiative spin-2 sector but LINEARLY and from a PRESCRIBED source, so the radiated energy is never taken back out "
"of the matter that emitted it. The three pieces -- relativistic quantum matter, the dynamical spin-2 field, and a "
"self-consistent exchange of energy between them -- each now exist separately; coupling all three in one evolution, so "
"that a source radiates and thereby decays, is the remaining integration problem, and it is a substantial build rather "
"than a further increment. Each barrier is met "
"individually, and two pairs now jointly; assembling them all into one consistent theory that also fixes the "
"Standard-Model content and the constants of nature is the work of fundamental physics itself, not of this toy. The value of Section 8 is to show that these barriers, usually treated as fatal to any "
"'space is a medium' program, are here concrete and — one at a time — surmountable. The quantum-mechanics result "
"of Section 8.6 carries its own explicit boundary: the wave, hbar, the Born statistics and de Broglie's relation "
"emerge, but the guidance of a particle by a separate pilot wave and the selection of a single definite outcome are "
"not derived — the measurement problem is left where it stands for every interpretation of quantum theory.")

# ===== 11 Conclusion =====
heading("11  Conclusion", 1)
body("Starting from space as an active medium, the project produced persistent particles, gave them charge "
"and spin as topological and Noether quantum numbers, showed those numbers robust on a self-organizing "
"substrate, and derived gravity from the medium — finding it real but short-ranged. Pursuing why, it "
"identified screening as a mass term and a conservation law as the escape, then realized, from one field "
"promoted step by step, a long-range Goldstone force, a mass-screened force, a gauge-screened (Meissner) "
"force, and finally a genuine 1/r² electromagnetism between quantized monopoles in three dimensions. The "
"model did not fail to make gravity Newtonian so much as explain, in its own terms, which forces are "
"long-range and why. It then turned the same measuring discipline on the question of whether it could be "
"fundamental, and found the classic obstacles — Lorentz invariance, chiral fermions, quantization, and even "
"spin-2 gravity — not fatal but, one at a time, met: an emergent Lorentz-invariant cone, a single chiral "
"fermion on a domain wall, and a relativistic quantum field theory. Pressed further, these separate victories "
"collapsed into one "
"principle. A field bolted on beside the medium brings a second light cone; a field made OF the medium inherits "
"the first. Carried to its end, the medium's own bond fluctuations, seen by its fermions, ARE the emergent "
"photon and the emergent graviton — so matter, light and gravity share a single cone not by tuning but by "
"construction. That graviton was then made dynamical: a luminal spin-2 wave whose Einstein-Hilbert kinetic term is "
"induced from the fermion loop rather than imposed. Turned on quantum mechanics "
"itself, the same discipline found most of it to be condensate mechanics: the Schrodinger wave and hbar as a "
"material property, the Born rule as a stochastic attractor, and de Broglie's lambda = h/p all emerge, with only "
"the guidance of a particle by a separate wave and the selection of a single outcome — the measurement problem's "
"hard core — left as a postulate. And, for the first time, the model made a prediction rather than a reproduction: "
"a Planck-suppressed, cross-species-universal, crystallographically-anisotropic quadratic Lorentz violation that "
"survives every present bound. Put against real numbers, however, that particular signature proved safe but out of "
"reach — a quadratic suppression sits some sixteen orders below current sensitivity — and the model's genuinely "
"reachable prediction turned out to be gravitational instead: a scale-dependent light-bending parameter below the "
"amplitude mode's Compton wavelength, in the sub-millimetre window short-range gravity experiments are entering now.")
body("The same discipline that produced those results also destroyed one, and then rebuilt it. Pressed on why gravity "
"cannot be shielded, the model gave a genuinely satisfying answer -- screening is neutralization, not loss, and a "
"topological charge is unneutralizable, so quantization does in the medium what the absence of negative mass does in "
"nature. But applying the same standard of proof to the FORCE, rather than to the field, brought the whole "
"gravitational programme down: two like curvature charges REPEL with a force that GROWS with separation, and the "
"tetrad graviton, for all its genuine long-range 1/r^2 field, exerts no long-range force at all. Every one of those "
"failures was an ELASTIC calculation, and naming that gave the answer. A mass is a force DIPOLE, so the elastic "
"sector cannot host gravity for any choice of moduli — while the principle the project had established for every "
"other force, that RANGE IS SET BY WHETHER A SYMMETRY PROTECTS THE MEDIATOR FROM A MASS TERM, had never once been "
"applied to gravity's own mediator. Applied at last, it identifies that mediator as the condensate's AMPLITUDE mode: "
"unprotected, therefore gapped, therefore Yukawa — which is exactly why every earlier measurement had found gravity "
"screened — and coupled monopolarly to positive-definite energy, so that its exchange between like charges is "
"universally ATTRACTIVE, with no sign inserted by hand. Its range is the inverse gap, measured; its potential is a "
"screening exponential times a 1/r Newtonian core, measured; and at criticality the exponential goes to unity and "
"Newton's law stands alone. Gravity, in this medium, is what the medium's own principle always said it would be.")
body("That scalar gravity was then carried the rest of the way to the tensor theory. The medium's curvature sector — "
"the graviton's true home — is CONFINING on its own, a string tension rather than a force; but the same fermion loop "
"that induces the photon induces an Einstein term that DECONFINES it into a Newtonian graviton, and the sign that "
"decides whether this works was MEASURED, positive, by holding the induced gravity against the model's own healthy "
"photon. In the physical four dimensions that graviton is a genuine dynamical spin-2 wave with two degenerate "
"helicity-2 polarizations, and its Einstein normalization — light bending by the famous factor of two — follows from "
"Weinberg's theorem as an INFRARED-emergent identity, in exactly the sense the model's Lorentz invariance is "
"emergent, since diffeomorphism invariance is not a symmetry of any lattice. The same condensate structure then "
"disarms the deepest quantitative disaster in physics: the vacuum's enormous zero-point energy does not gravitate, "
"because a self-sustained condensate gravitates its grand potential, which vanishes at equilibrium for ANY bare "
"value — the cosmological-constant fine-tuning dissolved, though not its observed residue derived. And the induction "
"mechanism that built the photon builds non-Abelian Yang-Mills as readily, for SU(2) and SU(3) alike, with chirality "
"riding alongside it consistently: the anomaly on a domain wall is a quantized integer, and the bulk supplies exactly "
"what each wall loses. The same confrontation with data that deflated the Lorentz prediction also retracted the "
"model's other great fine-tuning — the 1e122 approach to criticality, which had rested on the superseded scalar "
"reading of gravity — so that both of the program's 10^122 tunings have now fallen, each to a structural result "
"rather than a fitted parameter.")
body("What remains is honest and specific. General relativity is reached as an INFRARED FIXED POINT, not as a "
"lattice-exact law — the magnitude of Newton's constant is cutoff-dependent, and a direct lattice-exact measurement "
"of gamma = 1 is structurally unavailable because diffeomorphism invariance is emergent by construction. The observed "
"nonzero cosmological constant, the DERIVATION of the specific Standard-Model group and its representations and "
"constants (the group remains an input, only its Yang-Mills dynamics are induced), the residual measurement problem, "
"and quantitative rather than qualitative correspondences are, in the spirit of the project, measurements waiting to "
"be made. A model that can be made to say something false, and then be caught doing it, and then be made to say "
"something true, is the only kind worth building.")

# --- references ---
heading("References", 1)
refs=[
 "J. D. Eshelby, The continuum theory of lattice defects, Solid State Physics 3, 79 (1956). (Bitter–Crum; centres of dilatation.)",
 "J. Goldstone, Field theories with «superconductor» solutions, Nuovo Cimento 19, 154 (1961).",
 "J. M. Kosterlitz and D. J. Thouless, Ordering, metastability and phase transitions in two-dimensional systems, J. Phys. C 6, 1181 (1973).",
 "A. A. Abrikosov, On the magnetic properties of superconductors of the second group, Sov. Phys. JETP 5, 1174 (1957).",
 "P. W. Anderson, Plasmons, gauge invariance, and mass, Phys. Rev. 130, 439 (1963). (The Anderson–Higgs mechanism.)",
 "G. 't Hooft, Magnetic monopoles in unified gauge theories, Nucl. Phys. B 79, 276 (1974); A. M. Polyakov, JETP Lett. 20, 194 (1974).",
 "T. A. DeGrand and D. Toussaint, Topological excitations and Monte Carlo simulation of Abelian gauge theory, Phys. Rev. D 22, 2478 (1980).",
 "J. C. Maxwell, On the calculation of the equilibrium and stiffness of frames, Phil. Mag. 27, 294 (1864). (Central-force rigidity.)",
 "G. E. Volovik, The Universe in a Helium Droplet, Oxford (2003). (Emergent relativity, gauge fields and gravity near a Fermi point.)",
 "H. B. Nielsen and M. Ninomiya, Absence of neutrinos on a lattice, Nucl. Phys. B 185, 20 (1981). (Fermion doubling.)",
 "D. B. Kaplan, A method for simulating chiral fermions on the lattice, Phys. Lett. B 288, 342 (1992); C. G. Callan and J. A. Harvey, Nucl. Phys. B 250, 427 (1985). (Domain-wall chiral fermions.)",
 "M. Pretko and L. Radzihovsky, Fracton-elasticity duality, Phys. Rev. Lett. 120, 195301 (2018).",
 "S. Weinberg and E. Witten, Limits on massless particles, Phys. Lett. B 96, 59 (1980).",
 "A. D. Sakharov, Vacuum quantum fluctuations in curved space and the theory of gravitation, Dokl. Akad. Nauk SSSR 177, 70 (1967). (Induced gravity.)",
 "F. Guinea, M. I. Katsnelson and A. K. Geim, Energy gaps and a zero-field quantum Hall effect in graphene by strain engineering, Nature Phys. 6, 30 (2010). (Strain-induced emergent gauge fields.)",
 "K. G. Wilson, Confinement of quarks, Phys. Rev. D 10, 2445 (1974). (Lattice gauge theory; exact lattice gauge invariance and the Wilson plaquette action.)",
 "S. L. Adler, Einstein gravity as a symmetry-breaking effect in quantum field theory, Rev. Mod. Phys. 54, 729 (1982). (Induced gravity; the sign of the induced Newton constant.)",
 "S. Weinberg, Photons and gravitons in S-matrix theory: derivation of charge conservation and equality of gravitational and inertial mass, Phys. Rev. 135, B1049 (1964). (A massless spin-2 coupled to a conserved stress tensor is Einstein; gamma = 1.)",
 "R. Voss, Butler–Voss Condensate — project reference (CHEATSHEET) and source repository, 2026 (private).",
]
for i,rf in enumerate(refs,1):
    p=doc.add_paragraph(); p.paragraph_format.left_indent=Inches(0.3); p.paragraph_format.first_line_indent=Inches(-0.3)
    p.paragraph_format.space_after=Pt(4)
    nr=p.add_run(f"[{i}]  "); nr.bold=True; nr.font.size=Pt(9.5); nr.font.color.rgb=GREY
    p.add_run(rf).font.size=Pt(9.5)

apx=doc.add_paragraph(); ar=apx.add_run("Implementations (pure NumPy; private repository): simulation.py, "
"prototype_complex.py, h9_binding.py, prototype_mobile_nodes.py (H1–H10); integration_field_medium.py, "
"integration_phase2.py (§4); integration_phase3_variational.py, integration_phase3d.py, "
"density_response_3d_large.py (§5); screening_diagnosis.py, screening_gauss.py, screening_topocharge.py, "
"screening_gauged.py, screening_gauged_mobile.py (§6); screening_topocharge_3d.py, route_b_hedgehog.py, "
"route_c_monopole.py (§7); test_lorentz.py, test_lorentz_unified.py, test_lorentz_boost.py, test_dirac.py, "
"test_domain_wall.py, test_quantization.py, test_graviton.py, test_fracton_gravity.py, test_graviton_spin2.py, "
"test_cone_universality.py, test_cone_lock.py, test_induced_action.py, test_emergent_tetrad.py (§8.1-8.5); "
"test_emergent_qm.py, test_born_rule.py, test_double_solution.py (§8.6); test_graviton_dynamics.py, "
"test_induced_gravity.py (§8.7); test_lv_prediction.py (§8.8); test_shielding.py, test_tetrad_shielding.py, "
"test_disclination_force.py, test_tetrad_force.py (§8.9); test_critical_gravity.py (§8.10); test_two_gravities.py, "
"test_deconfinement.py (§8.11); test_induced_sign.py, test_spin2_dynamical.py, test_lattice_ward.py (§8.12); "
"test_cosmological_constant.py (§8.13); test_yang_mills.py (§8.14); test_anomaly_inflow.py (§8.15); "
"test_experimental_bounds.py (§8.16); "
"test_realtime_pump.py (§8.17); test_backreaction.py (§8.18); "
"test_newton_constant.py (§8.19); "
"test_gravitational_radiation.py (§8.20).")
ar.font.size=Pt(8.5); ar.font.color.rgb=GREY; ar.italic=True; apx.paragraph_format.space_before=Pt(12)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
doc.save(OUT)
if os.path.exists(OLD):
    os.remove(OLD)
print("SAVED:", OUT)
print("removed old screening-only doc:", not os.path.exists(OLD))
print("paragraphs:", len(doc.paragraphs), "tables:", len(doc.tables))
