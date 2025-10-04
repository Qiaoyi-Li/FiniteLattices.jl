using Documenter, FiniteLattices

makedocs(;
    modules = [FiniteLattices],
    sitename="FiniteLattices.jl",
    authors = "Qiaoyi Li",
    warnonly = [:missing_docs, :cross_references],
)

if haskey(ENV, "GITHUB_REF")
	@show ENV["GITHUB_REF"]
	branch = splitpath(ENV["GITHUB_REF"])[end]
	if branch == "main"
		devbranch = "main"
		devurl = "stable"
	elseif branch == "dev"
		devbranch = "dev"
		devurl = "dev"
	end
	deploydocs(
		repo = "github.com/Qiaoyi-Li/FiniteLattices.jl",
		devbranch = devbranch,
		devurl = devurl,
		push_preview = true,
	)
end