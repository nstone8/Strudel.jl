module Strudel
using Recluse

export writeconfig, scaffoldjob, scaffold

"""
```julia
writeconfig([filename])
```
Write a scaffold config file pre-filled with default values
to `filename` or `"config.jl"` if `filename is omitted.
"""
function writeconfig(filename="config.jl")
    configdict = Dict(
        :arraycenter => [0,0],
        :arraysize => [1,1],
        :arraypitch => [0,0],
        :calibratedfield => 1750,
        :laserpower => 100,
        :scanspeed => 100000,
        :stagevelocity => 100,
        :focusoffset => 100,
        :interfacepos => 1,
        :scafl => 6000,
        :scafw => 4000,

        :maxbeamlength => 120,
        :beamheight => 10,
        :beamwidth => 25,
        :beambottomz => 30,

        :overlapangle => (pi/6),
        :overlap => 10,
        :threadspeed => 10000,
        :segmentlength => 50,
        :gap => 15,
        :dhatch => 0.3,
        :dslice => 1,
        :dhammockhatch => 1,

        :postside => 100,
        :postheight => 80,
        :bumperfile => "bumper_data.gwl",
        :postkernelfile => "post_kernel_data.gwl",
        :kerneldir=>"kernels",
        :outputfile => "scaffoldjob.gwl",
        :backlashamount => 100,        
    )
    open(filename,"w") do io
        print(io,repr(configdict))
    end
end

"""
generate a `Vector` of the elements of a matrix with the elements reordered
to reflect a snake-style traversal rather than lexical
"""
function snakepattern(mat::Matrix)
    rows=map(1:size(mat)[2]) do rowindex
        this_row=mat[:,rowindex]
        isodd(rowindex) ? this_row : reverse(this_row)
    end
    vcat(rows...)
end

"""
```julia
scaffold(buffer[,centerposition],kwargs...)
```

Generate the files required to print a scaffold with geometry
defined by the provided `kwargs`. Append the top level job script
to the provided `buffer` (basic header entries should be added by
`scaffoldjob`). The scaffold will be centered on `centerposition`
which is `[0,0]` if not specified. All distances should be provided
in microns.

# Arguments
- bumperfile: .gwl file which will write the bumpers centered on `(0,0)`
- postkernelfile: .gwl file which will write one 'kernel' of posts, centered on `(0,0)`
- kerneldir: folder to create and fill with reusable units of the design
- calibratedfield: The diameter of the area which can be printed on without moving the stage
- threadspeed: scan speed for mesh
- interfacepos: how deep to start printing into the substrate
- scafl: scaffold length
- scafw: scaffold width
- maxbeamlength: maximum allowable length of connecting beams
- beamheight: height of connecting beams
- beamwidth: width of connecting beams
- beambottomz: z distance between the bottom of the beams and the bottom of the posts
- overlapangle: angle with which beams and beam segments should be chamfered
- overlap: overlap between neighboring connected geometry such as beam segments and beam/post connections
- segmentlength: maximum length of beam segments
- gap: gap that should be left between halves of beams during printing
- dhatch: hatching distance for solid objects
- dslice: slicing distance for solid objects
- dhammockhatch: hatching distance for mesh
- postside: side length of (square) posts
- postheight: post height
- backlashamount distance to move during backlash correction
- focusoffset: How far up should we move the objective in between scaffolds to ensure the autofocus works properly
"""
function scaffold(gwlbuf,centerposition=[0,0];calibratedfield,
                  interfacepos,
                  scafl,
                  scafw,
                  maxbeamlength,
                  beamheight,
                  beamwidth,
                  beambottomz,
                  overlapangle,
                  overlap,
                  threadspeed,
                  segmentlength,
                  gap,
                  dhatch,
                  dslice,
                  dhammockhatch,
                  postside,
                  postheight,
                  bumperfile,
                  postkernelfile,
                  kerneldir,
                  backlashamount,
                  focusoffset)
    
    hammockz = (beambottomz + beamheight/2)
    chamfer = [-overlapangle -overlapangle
	       overlapangle  overlapangle]
    #helper function to do backlash correction
    backlashstr =join([
        "MoveStageX $backlashamount",
        "MoveStageX $(-backlashamount)",
        "MoveStageY $backlashamount",
        "MoveStageY $(-backlashamount)"],"\n","\n")
    backlash(io) = println(io,backlashstr)

    bumperoffset = (scafw-postside)/2
    header = join([
        "GlobalGoToX $(centerposition[1])",
        "GlobalGoToY $(centerposition[2])",
        backlashstr,
        #move the objective up a bit before writing each scaffold to help autofocus
        "AddZDrivePosition $focusoffset",
        "FindInterfaceAt $interfacepos",
        "GlobalGoToX $(centerposition[1])",
        "GlobalGoToY $(centerposition[2]+bumperoffset)",
        backlashstr,
        "include $bumperfile",
        "AddZDrivePosition $(-1*postheight)",
        "GlobalGoToX $(centerposition[1])",
        "GlobalGoToY $(centerposition[2]-bumperoffset)",
        backlashstr,
        "include $bumperfile",
        "AddZDrivePosition $(-1*postheight)"
    ],"\n")
    
    println(gwlbuf,header) #pop our header into our buffer
    
    #calculate our grid dimensions
    nposty = ceil(Int, (scafw - 2*postside + maxbeamlength)/(postside+maxbeamlength))
    npostx = ceil(Int, (scafl + maxbeamlength)/(postside+maxbeamlength))

    #what are our 'actual' beam lengths?
    blx = (scafl - npostx*postside)/(npostx-1)
    bly = (scafw - nposty*postside - 2*postside)/(1+nposty)

    #define a helper function so we can get all multiples of our npost values
    function allmultiples(x)
	filter(1:x) do y
	    (x%y) == 0
	end
    end

    #what is the biggest 'kernel' of posts that 1) fits evenly into the desired space
    #and can fit inside a circle with diameter `calibratedfield`
    #we will assume each 'kernel' includes connecting bridges to the previous kernel
    #above and to the left

    allnumcombos = [(nx,ny) for nx in allmultiples(npostx), ny in allmultiples(nposty)]
    possiblenumcombos = filter(allnumcombos) do (nx,ny)
	sizex = nx*blx + nx*postside + overlap
	sizey = ny*bly + ny*postside + overlap
	sqrt(sizex^2 + sizey^2) < calibratedfield #does this fit
    end

    (_,maxcomboindex) = map(possiblenumcombos) do (nx,ny)
	nx*ny
    end |> findmax

    (nx,ny) = possiblenumcombos[maxcomboindex]

    #write out what we expect from our kernel file
    @info """
        The file $postkernelfile should contain a .gwl script which writes
        an array of posts with the following properties centered on (0,0).
        post side length: $postside
        post height: $postheight 
        array size: $nx × $ny
        array pitch in x direction : $(postside+blx)
        array pitch in y direction : $(postside+bly)
        """

    #size of kernels
    ksx = nx*blx + nx*postside + overlap
    ksy = ny*bly + ny*postside + overlap
    #calculate the coordinates of the center of each `kernel`, assuming the whole scaffold is centered at (0,0)
    #top left corner of the active area
    topleftx = -scafl/2
    toplefty = scafw/2 - postside
    firstkernelx = topleftx + ksx/2 - blx - overlap
    firstkernely = toplefty - ksy/2 + overlap
    kernelcenters = [(x,y) for x in firstkernelx:(ksx-overlap):(scafl/2), y in firstkernely:(-ksy+overlap):(-scafw/2 + postside + bly)]

    #translate these kernel centers to `centerposition`
    kernelcenters = map(kernelcenters) do kc
        (kc[1] + centerposition[1], kc[2] + centerposition[2])
    end
    #get the coordinates of each post relative to the top left of the kernel
    postoffsetstopleft = [(x,y) for x in range(start=overlap + blx + postside/2, step=blx+postside, length=nx),
		              y in range(start=-overlap - bly - postside/2, step=-bly-postside, length = ny)]
    #now relative to the center
    postoffsets = map(postoffsetstopleft) do (x,y)
	(x-(ksx/2),y+(ksy/2))
    end

    #write all the kernels
    for (kcx,kcy) in kernelcenters
	#we have to do different stuff if we are on the left, right or bottom edge
	leftkernel = kcx == kernelcenters[1,1][1]
	rightkernel = kcx == kernelcenters[end,end][1]
	bottomkernel = kcy == kernelcenters[end,end][2]
	#give each unique type of kernel a different name
	kernelname="kernel"
	if leftkernel
	    kernelname *= "left"	
	end
	
	if rightkernel
	    kernelname *= "right"
	end
	
	if bottomkernel
	    kernelname *= "bottom"
	end

	#move to kernel position (put this in main script so we can
	#reuse kernels
	println(gwlbuf,"GlobalGoToX $kcx")
	println(gwlbuf,"GlobalGoToY $kcy")
	backlash(gwlbuf)
        
	kernelfile=joinpath(kerneldir,"$kernelname.gwl")
        #include this file in our main job
	println(gwlbuf,"include $kernelfile")
        #don't regenerate a kernel we've already written
        isfile(kernelfile) ? continue : nothing
	open(kernelfile,"w") do kio
	    #the post array is not centered in this frame, add the offset
	    println(kio,"XOffset $((overlap+blx)/2)")
	    println(kio,"YOffset $((-overlap-bly)/2)")
	    #write the posts
	    println(kio,"include $(joinpath("..",postkernelfile))")
	    #undo the offset before we write all of our bridges, etc.
	    println(kio,"XOffset 0")
	    println(kio,"YOffset 0")
	    println(kio,"AddZDrivePosition $(-1*postheight)")
	    #write the posts
	    #all posts get bridges going up and to the left, except for posts along the left edge
	    #we will record all the hammocks that we need as we build the bridges and write them
	    #at the end of each kernel
	    hams = Vector{Hammock}()
	    #could throw the windows in with hams, but I don't want to
	    wins = Vector{Window}()
	    #offset between post centers and bridge centers
	    bridgeoffsets = [-1, 1]*(((postside-beamwidth)/2)-(beamheight/2)*tan(overlapangle))
	    for postcenter in postoffsets
		notleftedge = !leftkernel || (postcenter[1] != postoffsets[1,1][1])
		if notleftedge
		    #not on the left edge, write bridges to the left
		    bridgecenterx = postcenter[1] - (postside + blx)/2
		    bridgelengthx = blx + 2*overlap
		    #tan overlapangle term is to correct for the fact that the bottom of the beam is wider
		    #than beamwidth due to how chamfers are implemented
		    for bcy in postcenter[2] .+ bridgeoffsets
			b=Bridge(bridgelengthx,beamwidth,beamheight;
				 segmentlength, dslice, dhatch, overlap, overlapangle, gap, chamfer,
				 rotation=0, center=[bridgecenterx,bcy,hammockz])
			savegwl(kio,b)
		    end
		    #this pair of bridges gets a hammock
		    hamlength=sum(abs.(bridgeoffsets))-beamwidth+2*overlap
		    push!(hams,Hammock(hamlength,
				       bridgelengthx,hammockz,ceil(Int,hamlength/dhammockhatch),
				       rotation=pi/2,center=[bridgecenterx,postcenter[2]]))
		end
		#all posts get bridges directly above them
		bridgecentery = postcenter[2] + (bly + postside)/2
		bridgelengthy = bly + 2*overlap
		for bcx in postcenter[1] .+ bridgeoffsets
		    b=Bridge(bridgelengthy,beamwidth,beamheight;
			     segmentlength, dslice, dhatch, overlap, overlapangle, gap, chamfer,
			     rotation=pi/2, center=[bcx,bridgecentery,hammockz])
		    savegwl(kio,b)
		end
		#this pair of bridges gets a hammock
		hamwidth=sum(abs.(bridgeoffsets))-beamwidth+2*overlap
		push!(hams,Hammock(bridgelengthy,
				   hamwidth,hammockz,ceil(Int,bridgelengthy/dhammockhatch),
				   rotation=pi/2,center=[postcenter[1],bridgecentery]))
		
		#every post not on the left edge also gets a hammock diagonally to it's northwest
		notleftedge ? push!(hams,Hammock(
		    bly + 2*overlap, blx+2*overlap, hammockz,
		    ceil(Int,(bly+2*overlap)/dhammockhatch),
		    rotation=pi/2,center=[postcenter[1] - blx/2 - postside/2,
					  postcenter[2] + bly/2 + postside/2]
		)) : nothing
		#if we're on the bottom edge we need bridges directly below each post
		onbottomedge = bottomkernel && (postcenter[2] == postoffsets[end,end][2])
		if onbottomedge
		    #bridges directly below us
		    bridgecentery = postcenter[2] - (bly + postside)/2
		    bridgelengthy = bly + 2*overlap
		    for bcx in postcenter[1] .+ bridgeoffsets
			b=Bridge(bridgelengthy,beamwidth,beamheight;
				 segmentlength, dslice, dhatch, overlap, overlapangle, gap, chamfer,
				 rotation=pi/2, center=[bcx,bridgecentery,hammockz])
			savegwl(kio,b)
		    end
		    #this pair of bridges gets a hammock
		    hamwidth=sum(abs.(bridgeoffsets))-beamwidth+2*overlap
		    push!(hams,Hammock(bridgelengthy,
				       hamwidth,hammockz,ceil(Int,bridgelengthy/dhammockhatch),
				       rotation=pi/2,center=[postcenter[1],bridgecentery]))
		    
		    #every post not on the left edge also gets a hammock diagonally to it's southwest
		    notleftedge ? push!(hams,Hammock(
			bly + 2*overlap, blx+2*overlap, hammockz,
			ceil(Int,(bly+2*overlap)/dhammockhatch),
			rotation=pi/2,center=[postcenter[1] - blx/2 - postside/2,
					      postcenter[2] - bly/2 - postside/2]
		    )) : nothing
		end
		#if we're on the right or left edge we want Windows
		onrightedge = rightkernel && (postcenter[1] == postoffsets[end,end][1])
		onleftedge = !notleftedge
		if (onleftedge || onrightedge)
		    @assert onleftedge == !(onrightedge) "Sanity check"
		    winwidth = bly + 2*overlap
		    #get the z position of the bottom of the window
		    winz = hammockz + beamheight/2 - overlap
		    #need to divide by the cosine to compensate for the tilt
		    winheight = (postheight - winz)/cos(overlapangle)
		    winthreads = ceil(Int,winheight/dhammockhatch)
		    winx = onleftedge ? postcenter[1] - abs(bridgeoffsets[1]) :
			postcenter[1] + abs(bridgeoffsets[1])
		    winys=[postcenter[2]+(bly+postside)/2]
		    rot = onleftedge ? pi/2 : 3*pi/2
		    #also one underneat us if we're on the bottom edge
		    onbottomedge ? push!(winys,postcenter[2]-(bly+postside)/2) : nothing
		    for winy in winys
			push!(wins,Window(winwidth,winheight,winz,winthreads,
					  rotation=[rot,overlapangle],center=[winx,winy],
					  taper=overlapangle))
		    end
		end
	    end
	    #write the hammocks
	    println(kio,"ScanSpeed $threadspeed")
	    savegwl(kio,hams...)
	    #write the windows
	    savegwl(kio,wins...)
	end
    end
end


"""
```julia
scaffoldjob(;kwargs...)
scaffoldjob(configfile)
```

Generate the files required to print an array of scaffolds
with geometry defined by the provided `kwargs` or config file.
All arguments of `scaffold` are required, as well as:

# Additional arguments
- laserpower: laser power
- scanspeed: scan speed for solid objects
- stagevelocity: how fast should the stage move between kernels and scaffolds
- outputfile: where to write the toplevel job script
- arraycenter: where should the center of the scaffold array be in global coordinates, specified as a vector `[x,y]`
- arraysize: size of the array, specified as a vector `[m,n]`
- arraypitch: pitch of the array in x and y, specified as a vector `[m,n]`
"""
function scaffoldjob(;arraycenter,arraysize,arraypitch,laserpower,
                     scanspeed,kerneldir,outputfile,stagevelocity,
                     scaffoldargs...)
    #make a buffer to store our toplevel script
    gwlbuf=IOBuffer()
    #create our header
    header = join([
        "PowerScaling 1.0",
        "var \$solidLaserPower = $laserpower",
        "var \$solidScanSpeed = $scanspeed",
        "var \$baseLaserPower = \$solidLaserPower",
        "var \$baseScanSpeed = \$solidScanSpeed/2",
        "StageVelocity $stagevelocity"    
    ],"\n")
    #write the header to file
    println(gwlbuf,header)
    #create the kernel directory
    mkdir(kerneldir)

    #lay out our scaffold centers, if arraysize is [1,1], just return [[0,0]]
    scaf_uncentered=(prod(arraysize) > 1) ? [[x,y] for x in range(0,step=arraypitch[1],length=arraysize[1]),
                                               y in range(0,step=arraypitch[2],length=arraysize[2])] : [[0,0]]
    (max_xpos,max_ypos) = map(1:2) do i
        maximum(scaf_uncentered) do su
            su[i]
        end
    end

    #center the kernel
    scaf_centered = map(scaf_uncentered) do (sx,sy)
        [sx - max_xpos/2, sy - max_ypos/2]
    end

    #translate so the center of the array lies at `arraycenter`
    scaf_translated = map(scaf_centered) do sc
        sc + arraycenter
    end
    
    #write all of the scaffolds to gwlbuf and generate kernel files
    for st in snakepattern(scaf_translated)
        scaffold(gwlbuf,st;kerneldir,scaffoldargs...)
    end

    #write out toplevel script
    seek(gwlbuf,0)
    write(outputfile,read(gwlbuf))
end

function scaffoldjob(configfile)
    configdict = include(configfile)
    scaffoldjob(;configdict...)
end

end # module Strudel
