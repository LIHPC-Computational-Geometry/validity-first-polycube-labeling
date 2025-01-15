# https://github.com/charmbracelet/vhs

# Output file
Output automatic_polycube.gif

# Setup
Set Shell zsh # issue with bash currently https://github.com/charmbracelet/vhs/issues/419
Set FontSize 12
Set Width 600
Set Height 300
Set TypingSpeed 100ms
Set Padding 10
Set Margin 10
Set MarginFill "#6B50FF"
Set WindowBar Colorful
Set BorderRadius 10
Hide
Type 'cd ../build_Release'
Enter
Type 'clear'
Enter

# Recorded commands
Show
Type 'head ../data/B0/surface.obj'
Enter
Sleep 2s
Type './bin/automatic_polycube ../data/B0/surface.obj output=../data/B0/labeling.txt gui=false'
Enter
Sleep 2s
Type 'head ../data/B0/labeling.txt'
Enter
Sleep 2s

# Cleanup
Hide
Type 'rm ../data/B0/labeling.txt'
Enter