docker stop shiny
docker rm shiny
docker run --name shiny -d -p 3838:3838 -v ~/Desktop/virus_orf_finder_new:/srv/shiny-server/ yumapsc/vorffinder2.0
