from limb import load_james_limbs, load_kevin_limbs, load_heura_limbs2, ageInHours
from limb import getDist, reposition, spline_xy
from vedo import Spline, Line, procrustesAlignment, write, show, load
from vedo import printc, mag, colorMap, Axes, Text2D
from vedo.pyplot import histogram, plot
import numpy as np, os

################################################### procrustesAlignment
def align_proc(limbs, stages):

    averages = []

    #first iteration #################################
    filtered_limbs1 = []
    for istage in range(len(stages)):

        sel_limbs = []
        splines = []
        scores = []
        for lm in limbs:
            if stages[istage] not in lm.ageAsString: continue
            if lm.side != side[0] : continue

            sel_limbs.append(lm)
            spln = Spline(lm.datapoints, res=NPT)
            lm.Line = spln.c('black').alpha(0.4).lw(3)
            lm.ageIndex = istage
            splines.append(spln)

        procus = procrustesAlignment(splines, rigid=True)
        mean = procus.info['mean']
        procusLimbs = procus.unpack()

        CUT = 0.12 #######################

        for i, l in enumerate(procusLimbs):
            darr = mag(l.points()-mean)  # distance array
            score = np.sum(darr) / len(darr)/ l.averageSize()
            scores.append(score)
            limbs[i].LineReg = l
            if score > CUT:
                l.c('blue3')
            else:
                l.cmap('hot_r', darr)
                filtered_limbs1.append(sel_limbs[i])

        if viz>=2 and len(mean):
            lmean = Line(mean, lw=5, c='k')
            show([(lmean.clone().shift(0,0,1), *procusLimbs),
                 histogram(scores, c='g6', title=stages[istage],
                           xtitle=f"fit score (cut={CUT})")],
                 N=2, sharecam=0).screenshot(f"pics/procr1_E{stages[istage]}.png").close()

    #second iteration #################################
    filtered_limbs2 = []
    for stage in stages:

        sel_limbs = []
        splines = []
        scores = []
        for lm in filtered_limbs1:
            if stage not in lm.ageAsString:
                continue
            sel_limbs.append(lm)
            splines.append(lm.Line)

        procus = procrustesAlignment(splines, rigid=True)
        mean = procus.info['mean']
        procusLimbs = procus.unpack()

        CUT = 0.075 #######################

        for i, l in enumerate(procusLimbs):
            darr = mag(l.points()-mean)  # distance array
            score = np.sum(darr) / len(darr)/ l.averageSize()
            scores.append(score)
            if score > CUT:
                l.c('blue3')
            else:
                l.cmap('hot_r', darr)
                filtered_limbs2.append(sel_limbs[i])

        if viz>=2 and len(mean):
            lmean = Line(mean, lw=5, c='k')
            show([(lmean.clone().shift(0,0,1), *procusLimbs),
                 histogram(scores, c='g4', title=stage, xtitle=f"fit score (cut={CUT})")],
                 N=2, sharecam=0).screenshot(f"pics/procr2_{side}_E{stage}.png").close()


    #result #################################
    for stage in stages:

        sel_limbs = []
        splines = []
        scores = []
        for lm in filtered_limbs2:
            if stage not in lm.ageAsString:
                continue
            sel_limbs.append(lm)
            splines.append(lm.Line)

        procus = procrustesAlignment(splines, rigid=True)
        mean = procus.info['mean']
        procusLimbs = procus.unpack()

        for i, l in enumerate(procusLimbs):
            darr = mag(l.points() - mean)  # distance array
            score = np.sum(darr) / len(darr)/ l.averageSize()
            scores.append(score)
            l.cmap('hot_r', darr)

        if not len(mean):
            averages.append(None)
            continue

        lmean = Line(mean, lw=5, c='k')

        if viz>=2:
            show([(lmean.clone().shift(0,0,1), *procusLimbs),
                 histogram(scores, c='g2', title=stage, xtitle="fit score")],
                 N=2, sharecam=0).screenshot(f"pics/procr3_{side}_E{stage}.png").close()

        averages.append(lmean)  # APPEND RESULT

    return averages

#######################################
def align_icp(limbs, stages, averages=None):

    if averages is None:
        averages = align_proc(limbs, stages)

    assert len(averages) == len(stages)
    new_averages=[]

    for stage in stages:

        istage = stages.index(stage)
        refline = averages[istage]
        if not refline:
            printc("WARNING: average absent for stage E"+stage, c='r')
            continue

        # check if file exists
        npyfilename = f"stage_{side}_E{stage}_icp_{limbs[0].author}.npy"
        if os.path.isfile(npyfilename):
            printc("REUSING:", npyfilename, "(delete it to generate a new one)", invert=1)
            lines = load(npyfilename)
            assert lines[0].N() == NPT # if not must delete the npy file!

        else:

            r0 = refline.points()[0]
            r1 = refline.points()[-1]
            rm = (r1+r0)/2
            rv = (r1-r0)/mag(r1-r0)
            refavesize = refline.averageSize()

            lines = []
            scores = []
            printc("ICP fitting for stage E"+stage, side, limbs[0].author, c='y', end=' ')
            for lm in limbs:
                if stage not in lm.ageAsString: continue
                if lm.side != side[0] : continue

                printc('*', c='y', end='')

                p0 = lm.Line.points()[0]
                p1 = lm.Line.points()[-1]
                pm = (p1+p0)/2
                pv = (p1-p0)/mag(p1-p0)
                cr = np.cross(pv, rv)[2]
                if not cr:
                    angle = 0
                else:
                    angle = np.arccos(np.dot(rv,pv)) * 57.3 * cr/abs(cr)

                ln = lm.Line.clone()
                ln.origin(pm).shift(-pm).rotateZ(angle, locally=True).shift(+rm)
                ln.alignTo(refline, rigid=True, useCentroids=False)
                dist = (getDist(refline, ln) + getDist(ln, refline)) /2 /refavesize
                if np.isnan(dist):
                    continue
                scores.append(dist)
                lm.icp_score = dist

                if dist > CUT_ICP: continue    ######## CUT!

                col = colorMap(dist, "hot_r", 0, 1.)
                ln.c(col, 0.5).lw(2)
                lines.append(ln)

            lines = [Line(l) for l in lines] #bug?
            write(lines, npyfilename) #### DUMP FILE

            print()
            if viz>=1:
                show([(*lines, refline.clone().z(1), stage),
                       histogram(scores, xlim=(0,5), c='red4', bc='w',
                                 title=f"{side} Limb ICP score for E{stage}"),
                     ], N=2, sharecam=False, bg='bb', viewup='2d',
                ).screenshot(f"pics/stage_{side}_E{stage}_icp_{limbs[0].author}.png").close()

        printc("..working out the average for the current group of lines")

        reposition(lines) # rotate lines about horizontal with their base at (0,0)

        t, ts, allx, allxs, ally, allys = spline_xy(lines)

        figx  = plot(t, allx, 'g.', ma=0.1, title="x-spline")
        figx += plot(ts,allxs, "r-", lw=3)

        figy  = plot(t, ally, 'g.', ma=0.1, title="y-spline")
        figy += plot(ts,allys, "r-", lw=3)

        avepts = np.c_[allxs, allys, np.zeros(NPT)]
        spl = Spline(avepts, res=NPT).c('r5').lw(4)
        write(spl.clone(),
              f"E{stage}_{side}_{limbs[0].author}.vtk")  ### DUMP VTK line (final result)

        message = Text2D(f"E{stage} {side}", bg='k', c='blue5', s=1.5, pos="top-center")
        if viz>=1:
            show([ lines+[spl.clone().z(1), Axes(spl), message], figx, figy],
                 N=3, sharecam=0, viewup='2d', size=(2400,800),
                 ).screenshot(f"pics/stage_{side}_E{stage}_average.png").close()

        new_averages.append(spl.c(istage).flag(stage))

    return new_averages

###########################################################################################
###########################################################################################
viz  = 1
NPT = 200  # spline nr of points
CUT_ICP = 0.5

############# process heura data
# side = "LEFT"
# averages = align_icp( * load_heura_limbs2("data/limbs_db.npy") )
# show(averages, f"Heura's dataset\n{side} limbs",
#      viewup='2d', axes=1).screenshot(f"pics/pic_{side}_heura.png").close()
# side = "RIGHT"
# averages = align_icp( * load_heura_limbs2("data/limbs_db.npy") )
# show(averages, f"Heura's dataset\n{side} limbs",
#      viewup='2d', axes=1).screenshot(f"pics/pic_{side}_heura.png").close()

############# process james data
# side = "LEFT"
# averages = align_icp( * load_james_limbs("data/limbs_db.npy") )
# show(averages, f"James' dataset\n{side} limbs",
#      viewup='2d', axes=1).screenshot(f"pics/pic_{side}_james.png").close()
side = "RIGHT"
averages = align_icp( * load_james_limbs("data/limbs_db.npy") )
show(averages, f"James' dataset\n{side} limbs",
      viewup='2d', axes=1).screenshot(f"pics/pic_{side}_james.png").close()

############# process kevin data
# side = "LEFT"
# averages = align_icp( * load_kevin_limbs("data/limbs_db.npy") )
# show(averages, f"Kevin's dataset\n{side} limbs",
#      viewup='2d', axes=1).screenshot(f"pics/pic_{side}_kevin.png").close()
# side = "RIGHT"
# averages = align_icp( * load_kevin_limbs("data/limbs_db.npy") )
# show(averages, f"Kevin's dataset\n{side} limbs",
#       viewup='2d', axes=1).screenshot(f"pics/pic_{side}_kevin.png").close()



















