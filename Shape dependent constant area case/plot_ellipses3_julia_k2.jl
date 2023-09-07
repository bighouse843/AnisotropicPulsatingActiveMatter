using Base.Threads
using Distributed
using FFMPEG
using DelimitedFiles
using Plots
using LaTeXStrings

#In this code we visualize the results of the simulations performed using the C++ code. This code is meant to be run on multiple threads. 

addprocs(8) #8 is the number of threads and it cand be changed
@everywhere begin
    using PyCall
    using PyPlot
    using ProgressMeter

    np = pyimport("numpy")
    patch = pyimport("matplotlib.patches")
    cm = pyimport("matplotlib.cm")
    colors = pyimport("matplotlib.colors")
    half_cmap = cm["hsv"]
    new_colors = np.vstack((half_cmap(np.linspace(0, 1, 128)), half_cmap(np.linspace(0, 1, 128))))
    cmap = cm["hsv"]
    newcmap = colors.ListedColormap(new_colors, name="OrangeBlue")
    cmap1 = cm.get_cmap(newcmap)

    d = parse.(Float64, readlines("test"))
    Lx       = d[1]
    Ly       = d[2]
    ρ     = d[3]
    T        = d[4]
    omega    = d[8]
    t_total  = d[10]
    t_record = d[12]
    radius   = d[14]
    amp      = d[16]
    λ      = d[18]


    function read_n_lines(filename)
        open(filename) do file
            arr = readlines(file)
        end
        arr = split.(arr, "\t")
        arr2 = Array{Any}(undef, length(arr))
        for i in range(1,length(arr))
            arr2[i] = parse.(Float64, arr[i])
        end
        return arr2
    end
    #=

    function add_ellips(pos)
        x, y, ϕ, θ =  pos
        w = 2*radius*(1+λ*sin(ϕ))/(1+λ)
        ang = (θ+π/2)* 57.2958
        #println(typeof(y))
        c = ϕ%(2*π)/(2*π)
        ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
        axs.add_artist(ells)
        ells.set_clip_box(axs.bbox)
        ells.set_facecolor(cmap(c))
        ells.set_alpha(α) #trasparenza
        return
    end
    =#
    
    α = 0.8
    a = sqrt(2/(sqrt(3)*ρ))
    Np = Int(2*round(Lx/a)*round(Ly/(sqrt(3)*a)))
    nb_record = Int(min(1e3,t_total/t_record)) #20
    radius_2 = 2*radius
    π2 = 2*π
    π_2 = π/2

    
    #println(lista_file)
    function crea_singolo_frame_fase()
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        fig.suptitle("Phase", fontsize=30)
        axs.set_title("ρ="*string(ρ)*", ε="*string(amp)*", λ="*string(λ)*", ω="*string(omega)*", Lx="*string(Lx)*", Ly="*string(Ly)*", t="*string(omega*f*t_record/π2)[1:5])
        
        
        data = read_n_lines("test_snap")
    
        @inbounds for i in range(1,Np)
            x, y, ϕ, θ =  data[i]
            sin_lam = λ*sin(ϕ)
            w = radius_2*exp(sin_lam)
            h = radius_2*exp(-sin_lam)
            ang = (θ+π_2)* 57.2958
            #println(typeof(y))
            c1 = mod(ϕ,π2)/(π2)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap1(c1))
            ells.set_alpha(α) #trasparenza

        end
        fig.show()
        fig.savefig(string(0))
        plt.close()      
        return
    end
    function crea_frame_fase(f)
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        fig.suptitle("Phase", fontsize=30)
        axs.set_title("ρ="*string(ρ)*", ε="*string(amp)*", λ="*string(λ)*", ω="*string(omega)*", Lx="*string(Lx)*", Ly="*string(Ly)*", t="*string(omega*f*t_record/π2)[1:5])
        
        data = read_n_lines("test_snap_"*string(f))
    
        @inbounds for i in range(1,Np)
            x, y, ϕ, θ =  data[i]
            sin_lam = λ*sin(ϕ)
            w = radius_2*exp(sin_lam)
            h = radius_2*exp(-sin_lam)
            ang = (θ+π_2)* 57.2958
            #println(typeof(y))
            c1 = mod(ϕ,π)/(π)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap1(c1))
            ells.set_alpha(α) #trasparenza

        end
        fig.savefig(string(f))
        plt.close()      
        return
    end

 


    function crea_frame_orient(f)
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        fig.suptitle("Orientation", fontsize=30)
        axs.set_title("ρ="*string(ρ)*", ε="*string(amp)*", λ="*string(λ)*", ω="*string(omega)*", Lx="*string(Lx)*", Ly="*string(Ly)*", t="*string(omega*f*t_record/π2)[1:5])
        
        data = read_n_lines("test_snap_"*string(f))
    
        @inbounds for i in range(1,Np)
            x, y, ϕ, θ =  data[i]
            sin_lam = λ*sin(ϕ)
            w = radius_2*exp(sin_lam)
            h = radius_2*exp(-sin_lam)
            ang = (θ+π_2)* 57.2958
            
            θ = θ + π_2*ceil(ϕ/π)
            c2 = mod(θ,π)/(π)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap(c2))
            ells.set_alpha(α) #trasparenza
        end
        fig.savefig(string(f))
        plt.close()      
        return
    end

    function crea_frame_vel(f)
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        axs.set_title("Orientation of velocity", fontsize=30)
        
        
        data = read_n_lines("test_snap_"*string(f))
        data0 = read_n_lines("test_snap_"*string(f-10))
    
        @inbounds for i in range(1,Np)
            x, y, ϕ, θ =  data[i]
            x0, y0, ϕ0, θ0 =  data0[i]
            sin_lam = λ*sin(ϕ)
            w = radius_2*exp(sin_lam)
            h = radius_2*exp(-sin_lam)
            ang = (θ+π_2)* 57.2958
            #println(typeof(y))
            temp = atan(y-y0, x-x0)
            c1 = mod(temp,π2)/(π2)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap(c1))
            ells.set_alpha(α) #trasparenza



        end
        fig.savefig(string(f-10))
        plt.close()      
        return
    end

    function crea_frame_orient_change(f)
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        axs.set_title("Change of orientation", fontsize=30)
        
        data0 = read_n_lines("test_snap_"*string(f-10))
        data = read_n_lines("test_snap_"*string(f))

    
        @inbounds for i in range(1,Np)
            x, y, ϕ, θ =  data[i]
            x0, y0, ϕ0, θ0 =  data0[i]
            sin_lam = λ*sin(ϕ)
            w = radius_2*exp(sin_lam)
            h = radius_2*exp(-sin_lam)
            ang = (θ+π_2)* 57.2958
            #println(typeof(y))
            c1 = mod(θ-θ0,π)/(π)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap(c1))
            ells.set_alpha(α) #trasparenza



        end
        fig.savefig(string(f-10))
        plt.close()      
        return
    end

    function crea_frame_f11(f)
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        axs.set_title("f11", fontsize=30)
        
        
        data = read_n_lines("test_snap_"*string(f))
    
        @inbounds for i in range(1,Np)
            x, y, ϕ, θ =  data[i]
            sin_lam = λ*sin(ϕ)
            w = radius_2*exp(sin_lam)
            h = radius_2*exp(-sin_lam)
            ang = (θ+π_2)* 57.2958
            c2 = mod(2*θ+ϕ,π2)/(π2)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap(c2))
            ells.set_alpha(α) #trasparenza
        end
        fig.savefig(string(f))
        plt.close()      
        return
    end

    function crea_frame_f1_1(f)
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        axs.set_title("f1,-1", fontsize=30)
        
        
        data = read_n_lines("test_snap_"*string(f))
    
        @inbounds for i in range(1,Np)
            x, y, ϕ, θ =  data[i]
            sin_lam = λ*sin(ϕ)
            w = radius_2*exp(sin_lam)
            h = radius_2*exp(-sin_lam)
            ang = (θ+π_2)* 57.2958
            c2 = mod(2*θ-ϕ,π2)/(π2)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap(c2))
            ells.set_alpha(α) #trasparenza
        end
        fig.savefig(string(f))
        plt.close()      
        return
    end

    function crea_frame_f20(f)
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        axs.set_title("f20", fontsize=30)
        
        
        data = read_n_lines("test_snap_"*string(f))
    
        @inbounds for i in range(1,Np)
            x, y, ϕ, θ =  data[i]
            sin_lam = λ*sin(ϕ)
            w = radius_2*exp(sin_lam)
            h = radius_2*exp(-sin_lam)
            ang = (θ+π_2)* 57.2958
            c2 = mod(2*ϕ,π2)/(π2)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap(c2))
            ells.set_alpha(α) #trasparenza
        end
        fig.savefig(string(f))
        plt.close()      
        return
    end

    function crea_frame_f02(f)
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        axs.set_title("f02", fontsize=30)
        
        
        data = read_n_lines("test_snap_"*string(f))
    
        @inbounds for i in range(1,Np)
            x, y, ϕ, θ =  data[i]
            sin_lam = λ*sin(ϕ)
            w = radius_2*exp(sin_lam)
            h = radius_2*exp(-sin_lam)
            ang = (θ+π_2)* 57.2958
            c2 = mod(4*θ,π2)/(π2)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap(c2))
            ells.set_alpha(α) #trasparenza
        end
        fig.savefig(string(f))
        plt.close()      
        return
    end

end



#foreach(rm, filter(endswith(".png"), readdir())) #cancel old frames
#crea_singolo_frame_fase()

lista_file = filter(x->occursin("test_snap_",x), readdir())
nb_snaps = length(lista_file)

@showprogress pmap(crea_frame_orient, range(1,nb_snaps-1)) #here the functions presented above can be called. 
#crea_frame_fase creates the visualization of the phase of the particles
#crea_frame_orient creates the visualization of the orientation of the particles
#crea_frame_f11 creates the visualization of the phase plus the orientation

video_name = "rho_"*string(ρ)*"_eps_"*string(amp)*"phase"
#@ffmpeg_env run(`$ffmpeg -framerate 25 -i %d.png -c:v libx264 -profile:v high -crf 20 -y -pix_fmt yuv420p $video_name.mp4`)
#=
foreach(rm, filter(endswith(".png"), readdir())) #cancel old frames
@showprogress pmap(crea_frame_orient, range(1,nb_snaps-1))
video_name = "rho_"*string(ρ)*"_eps_"*string(amp)*"orientation"
@ffmpeg_env run(`$ffmpeg -framerate 25 -i %d.png -c:v libx264 -profile:v high -crf 20 -y -pix_fmt yuv420p $video_name.mp4`)

foreach(rm, filter(endswith(".png"), readdir())) #cancel old frames
@showprogress pmap(crea_frame_vel, range(11,nb_snaps-1))
video_name = "rho_"*string(ρ)*"_eps_"*string(amp)*"orientation of velocity"
@ffmpeg_env run(`$ffmpeg -framerate 25 -i %d.png -c:v libx264 -profile:v high -crf 20 -y -pix_fmt yuv420p $video_name.mp4`)

foreach(rm, filter(endswith(".png"), readdir())) #cancel old frames
@showprogress pmap(crea_frame_orient_change, range(11,nb_snaps-1))
video_name = "rho_"*string(ρ)*"_eps_"*string(amp)*"change of orientation"
@ffmpeg_env run(`$ffmpeg -framerate 25 -i %d.png -c:v libx264 -profile:v high -crf 20 -y -pix_fmt yuv420p $video_name.mp4`)


=#

#Per creare il video, usa il seguente comando da terminale
#ffmpeg -framerate 25 -i %d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4


    

#plt.show()




