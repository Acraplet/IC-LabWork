for f in ../data/data_powermeter_3*.csv; do
    python analyse_powermeter.py "$f"
done
