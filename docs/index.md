# Welcome to MkDocs

For full documentation visit [mkdocs.org](http://mkdocs.org).

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs help` - Print this help message.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.


## work with Docile.jl and Lexicon.jl

using Lexicon
index = save("docs/api/Lexicon.md", Lexicon);
save("docs/api/index.md", Index([index]); md_subheader = :category);
run(`mkdocs gh-deploy --clean`)