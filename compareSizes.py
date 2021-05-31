from limb import load_james_limbs, load_kevin_limbs, load_heura_limbs2, distPD
from vedo import Text3D, show
from vedo.pyplot import plot, fit

xlim = (220, 400) # plot limits
ylim = (  0, 950)


#####################################################################################
heu_limbs, _ = load_heura_limbs2("data/limbs_db.npy")
print("#heuras limbs", len(heu_limbs))

avsiz, times = [], []
for l in heu_limbs:
    avsiz.append(distPD(l.datapoints, isHeura=True)/1.2)
    times.append(l.age-10)

plt = plot(times, avsiz, '*', ms=0.8, mc='pink6', ma=0.8,
           title="Proximo-Distal length vs time",
           xlim=xlim, ylim=ylim,
           xtitle="time /h", ytitle="length /\mum",
)

txth = Text3D("Heura", [230,190], c='pink6', s=20).rotateZ(30, locally=True)
pfit = fit([times, avsiz], deg=2, niter=500, nstd=5, c='pink6')
plt += [txth, pfit, pfit.errorBand, *pfit.errorLines] # add these objects to Plot

#####################################################################################
james_limbs, _ = load_james_limbs("data/limbs_db.npy")
print("#james limbs", len(james_limbs))

avsiz, times = [], []
for l in james_limbs:
    avsiz.append(distPD(l.datapoints))
    times.append(l.age)

plt += plot(times, avsiz, 'a', xlim=xlim, ylim=ylim, ms=0.7, ma=0.3, mc='b4')

txtj = Text3D("James", [260,400], c='b3', s=20).rotateZ(30, locally=True)
pfit = fit([times, avsiz], deg=3, niter=500, nstd=5, c='b4')
plt += [txtj, pfit, pfit.errorBand, *pfit.errorLines] # add these objects to Plot

#####################################################################################
kevin_limbs, _ = load_kevin_limbs("data/limbs_db.npy")
print("#kevin limbs", len(kevin_limbs))

avsiz, times = [], []
for l in kevin_limbs:
    avsiz.append(distPD(l.datapoints)/1.7)
    if 290<l.age<300: l.age += 1.5
    times.append(l.age)

plt += plot(times, avsiz, '.', xlim=xlim, ylim=ylim, ms=1.1, ma=0.4, mc='g3')

txtk = Text3D("Kevin", [340,700], c='g3', s=20).rotateZ(30, locally=True)
pfit = fit([times, avsiz], deg=4, niter=500, nstd=5, c='g3')
plt += [txtk, pfit, pfit.errorBand, *pfit.errorLines] # add these objects to Plot

show(plt, size=(1800,1285), zoom=1.5, viewup="2d").screenshot('pics/sizes_time.png').close()


