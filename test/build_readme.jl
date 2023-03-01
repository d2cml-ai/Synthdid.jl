using Weave, Highlights


file = "README.jmd"

# @time weave(hdm_jl; doctype = "md2pdf", highlight_theme = Highlights.Themes.MonokaiMiniTheme, pandoc_options = ["--toc", "--toc-depth= 3", "--number-sections", "--self-contained"])
@time weave(file; doctype="github", pandoc_options=["--toc", "--toc-depth= 1", "--number-sections", "--self-contained"], highlight_theme=Highlights.Themes.TangoTheme)
# @time weave(hdm_jl; doctype="md2pdf", pandoc_options=["--toc", "--toc-depth= 1", "--number-sections", "--self-contained"])


