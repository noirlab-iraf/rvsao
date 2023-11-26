# File saotdc/pkg/untools/relearn.cl
# By Doug Mink, Center for Astrophysics
# November 17, 1992

procedure relearn (taskname)

char taskname

begin

    if (access ("./relearn.tmp"))
	delete ("./relearn.tmp")
    dpar (taskname,>"./relearn.tmp")
    unlearn (taskname)
    cl (<"./relearn.tmp")
    delete ("./relearn.tmp")

end
