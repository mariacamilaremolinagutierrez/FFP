all: run.py
	amuse run.py

clean:
	rm status.txt
	rm results.txt
	rm *.png
