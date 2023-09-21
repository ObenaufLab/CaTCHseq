import numpy as np
import logomaker as lm

 def generateLogo(self) -> 'logomaker.src.Logo.Logo':
        counts = lm.transform_matrix(pd.DataFrame(self._counts), from_type = "counts", to_type = "information")
        logo_obj = lm.Logo(counts, color_scheme = "classic")
        return logo_obj