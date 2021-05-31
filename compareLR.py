from vedo import Mesh, merge, show

def plot_lr(author, stage):
    l1 = Mesh(f"{path}/E{stage}_LEFT_{author}.vtk").c('r5').lw(4)
    l2 = Mesh(f"{path}/E{stage}_RIGHT_{author}.vtk").c('b5').lw(4)
    l2.alignTo(l1, rigid=1)
    ref=None
    # if author != "heura":
    #     ref= Mesh(f"../References/E{stage.replace('.',';')}_Reference.vtk")
    #     ref.mirror('y').c('k2').lw(4)
    #     s0 = ref.averageSize()
    #     s1 = l1.averageSize()
    #     s2 = l2.averageSize()
    #     s = (s1+s2)/2
    #     ref.scale(s/s0)
    #     ref.alignTo(merge(l1,l2), rigid=1)

    message = f"E{stage} by {author} \nred=LEFT\nblue=RIGHT\nblack=S.Sys."
    return show(l1,l2,ref, message, axes=1).screenshot(f"pics/LR_E{stage}_{author}.png")

path = "output"

plot_lr("heura", "10.00").close()
plot_lr("heura", "10.03").close()
plot_lr("heura", "10.06").close()
plot_lr("heura", "10.11").close()
plot_lr("heura", "10.14").close()

plot_lr("james", "10.21").close()
plot_lr("james", "11.09").close()
plot_lr("james", "11.21").close()
plot_lr("james", "12.09").close()
plot_lr("james", "12.21").close()

plot_lr("kevin", "12.09").close()
plot_lr("kevin", "13.09").close()
plot_lr("kevin", "14.09").close()
plot_lr("kevin", "14.21").close()
plot_lr("kevin", "15.09").close()
plot_lr("kevin", "15.21").close()
