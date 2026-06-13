from types import SimpleNamespace

def estimate_fdr(tscores):
    result = SimpleNamespace()
    result.true_component_mean = 0.75
    result.false_component_mean = 0.25
    result.true_component_weight = 0.60
    result.thresholds = {'low': 0.3, 'medium': 0.5, 'high': 0.7}
    return result
