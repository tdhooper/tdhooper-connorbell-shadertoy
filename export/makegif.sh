
#!/bin/sh

# From http://blog.pkh.me/p/21-high-quality-gif-with-ffmpeg.html

# Usage:
# ./makegif.sh "img%02d.png" output.gif

palette="/tmp/palette.png"

filters="scale=800:-1:flags=lanczos"

ffmpeg -v warning -start_number 1 -i $1 -vf "$filters,palettegen=stats_mode=diff:reserve_transparent=off:max_colors=256" -y $palette
ffmpeg -v warning -start_number 1 -i $1 -i $palette -lavfi "$filters [x]; [x][1:v] paletteuse" -y $2
