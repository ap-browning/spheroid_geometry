gr()
default()
default(
    fontfamily="Arial",
    tick_direction=:out,
    guidefontsize=9,
    annotationfontfamily="Arial",
    annotationfontsize=10,
    annotationhalign=:left,
    box=:on,
    msw=0.0
)

alphabet = "abcdefghijklmnopqrstuvwxyz"

function add_plot_labels!(plt;offset=0)
    n = length(plt.subplots)
    for i = 1:n
        plot!(plt,subplot=i,title="($(alphabet[i+offset]))")
    end
    plot!(
        titlelocation = :left,
        titlefontsize = 10,
        titlefontfamily = "Arial"
    )
end