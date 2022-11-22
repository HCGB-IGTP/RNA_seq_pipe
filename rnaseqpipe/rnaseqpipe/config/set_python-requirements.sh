pip list --format freeze | sed 's/\=\=/,/' > python/python_requirements.csv
awk -F',' '{print $1}' python/python_requirements.csv > python/python_requirement_summary.txt
