import resdk
from resdk.tables import QCTables

app = resdk.Resolwe(url="https://app.genialis.com")

# c = app.collection.get("fusions-seraseq-fusion-v3")
# c = app.collection.get("ddr-lee-et-al-2020-gse126765")
c = app.collection.get("qctables-test-collection")

qt = QCTables(c)
qt.clear_cache()
qt = QCTables(c)

# ----
print("=" * 100)
print(qt.general)
print(qt.qorts)
print(qt.qc)
