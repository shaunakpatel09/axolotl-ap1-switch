import warnings; warnings.filterwarnings('ignore')
import gsa
for a, b in [(10, 100), (100, 200), (200, 300), (300, 400)]:
    gsa.run(a, b)
print('ALL DONE', flush=True)
