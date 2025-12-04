for f in ../data/OPM/OPM_datasets_04Dec25/data_powermeter_*.csv; do
    python analyse_powermeter.py "$f"
done
