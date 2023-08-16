using Base.Threads
using Distributed
using FFMPEG
using DelimitedFiles
using Plots
using LaTeXStrings

addprocs(8)
@everywhere begin
    using PyCall
    using PyPlot
    using ProgressMeter

    patch = pyimport("matplotlib.patches")
    cm = pyimport("matplotlib.cm")
    cmap = cm["hsv"]

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

    function read_n_lines(filename, N)
        arr = Array{Any}(undef, N)
        open(filename) do file
            arr = readlines(file)
        end
        arr = split.(arr, "\t")
        arr2 = Array{Any}(undef, N)
        for i in range(1,N)
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
    w = 2*radius 
    α = 0.8
    a = sqrt(2/(sqrt(3)*ρ))
    Np = Int(2*round(Lx/a)*round(Ly/(sqrt(3)*a)))
    nb_record = Int(min(1e3,t_total/t_record)) #20
    radius_2 = 2*radius
    λ_1 = 1-λ
    π2 = 2*π
    π_2 = π/2

    
    #println(lista_file)

    function crea_frame_fase(f)
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        #axs.clear()
        #println("sono vivo1")
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        fig.suptitle("Phase", fontsize=30)
        axs.set_title("ρ="*string(ρ)*", ε="*string(amp)*", λ="*string(λ)*", ω="*string(omega)*", Lx="*string(Lx)*", Ly="*string(Ly)*", t="*string(omega*f*t_record/π2)[1:5])
        
        data = read_n_lines("test_snap_"*string(f), Np)
        #end
        #pmap(add_ellips, data)
    
        @inbounds for i in range(1,Np)
            #add_ellips(data[i])
            x, y, ϕ, θ =  data[i]
            h = radius_2*(1+λ*sin(ϕ))/λ_1
            ang = (θ+π_2)* 57.2958
            #println(typeof(y))
            c = mod(ϕ,π2)/(π2)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap(c))
            ells.set_alpha(α) #trasparenza
        end
        fig.savefig(string(f))
        plt.close()      
        return
    end

    function crea_frame_orient(f)
        fig, axs = plt.subplots(1, 1,figsize=(6,6))
        #axs.clear()
        #println("sono vivo1")
        axs.axis("equal")
        axs.set_xlim(0,Lx)
        axs.set_ylim(0,Ly)
        axs.set_xticks([])
        axs.set_yticks([])
        fig.suptitle("Orientation", fontsize=30)
        axs.set_title("ρ="*string(ρ)*", ε="*string(amp)*", λ="*string(λ)*", ω="*string(omega)*", Lx="*string(Lx)*", Ly="*string(Ly)*", t="*string(omega*f*t_record/π2)[1:5])
        
        data = read_n_lines("test_snap_"*string(f), Np)
        #end
        #pmap(add_ellips, data)
    
        @inbounds for i in range(1,Np)
            #add_ellips(data[i])
            x, y, ϕ, θ =  data[i]
            h = radius_2*(1+λ*sin(ϕ))/λ_1
            ang = (θ+π_2)* 57.2958
            #println(typeof(y))
            c = mod(θ,π)/(π)
            ells = patch.Ellipse(xy=(x,y), width=w, height=h, angle= ang)
            axs.add_artist(ells)
            ells.set_clip_box(axs.bbox)
            ells.set_facecolor(cmap(c))
            ells.set_alpha(α) #trasparenza
        end
        fig.savefig(string(f))
        plt.close()      
        return
    end

end
foreach(rm, filter(endswith(".png"), readdir())) #cancel old frames
lista_file = filter(x->occursin("test_snap_",x), readdir())
nb_snaps = length(lista_file)
@showprogress pmap(crea_frame_fase, range(1,nb_snaps-1))
video_name = "rho_"*string(ρ)*"_eps_"*string(amp)
@ffmpeg_env run(`$ffmpeg -framerate 25 -i %d.png -c:v libx264 -profile:v high -crf 20 -y -pix_fmt yuv420p $video_name.mp4`)



#Per creare il video, usa il seguente comando da terminale
#ffmpeg -framerate 25 -i %d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4


    

#plt.show()




