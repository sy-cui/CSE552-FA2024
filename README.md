# CSE552-FA2024
CSE 552: Nonlinear Finite Element Methods

## File structure of the Latex source code

```
📦doc
 ┣ 📂homework
 ┃ ┣ 📂hw1
 ┃ ┃ ┣ 📜hw1.tex
 ┃ ┃ ┣ 📜image_file_1.pdf
 ┃ ┃ ┣ ...
 ┃ ┣ 📂hw2
 ┃ ┣ 📂hw3
 ┃ ┗ ...
 ┣ 📂midterm
 ┃ ┣ 📜midterm.tex
 ┃ ┣ 📜midterm_image_file_1.pdf
 ┃ ┗ ...
 ┣ 📜cse552.tex
 ┣ 📜_settings.tex
 ┣ 📜cse552.pdf
 ┗ (other LATEX auxiliary files)
```

The main source code for each homework / exam are located in the `.tex` file in their respective directories. For instance, `/doc/midterm/midterm.tex` contains the source code for the midterm exam. 

### `cse552.tex`
Main driver file for generating the document. `PDFLATEX` is used for compilation. To select which homework / project to generate, `input` the source file from the respective directory and comment out the rest, e.g. 
```latex
% Generate pdf for the midterm exam
% \input{homework/hw1/hw1}
% \input{homework/hw2/hw2}
% ...

\input{midterm/midterm.tex}
```

### `_settings.tex`
Include preambles for custom settings / commands. For example, all `\newcommand` are defined here. 
