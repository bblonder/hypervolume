hypervolume_save_animated_gif <- function(image.size=400, axis=c(0,0,1),rpm=4,duration=15,fps=10,file.name='movie',directory.output='.',...)
{
  td = tempdir()
  tf = basename(tempfile(tmpdir=td))
  
  rgl::par3d(windowRect=c(0,0,image.size,image.size))
  rgl::movie3d(spin3d(axis=axis,rpm=rpm),duration=duration,fps=fps,movie=tf,dir=td,...)
  
  if(!file.exists(directory.output))
  {
    dir.create(directory.output)
  }
  file.rename(sprintf("%s.gif",file.path(td,tf)),file.path(directory.output,sprintf("%s.gif",file.name)))
}