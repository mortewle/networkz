
# hindre U-sving rett etter et svingforbud?
# vanskelig å få til uten å gjøre ting veldig komplisert og sikkert treigt

def u_sving():
    if u_sving_lov is False:
        from math import atan2, degrees
        def angle_between(p1, p2, p3):
            x1, y1 = p1.xy
            x2, y2 = p2.xy
            x3, y3 = p3.xy
            deg1 = (360 + degrees(atan2(x1[0] - x2[0], y1[0] - y2[0]))) % 360
            deg2 = (360 + degrees(atan2(x3[0] - x2[0], y3[0] - y2[0]))) % 360
            return deg2 - deg1 if deg1 <= deg2 else 360 - (deg1 - deg2)
        for p1, p2 in zip(dobbellenker["source_wkt"], dobbellenker["target_wkt"]):
            fra = veger_edges[veger_edges["source_wkt"]==p2]
            for idx, p3 in zip(fra.index, fra["target_wkt"]):
                sving_grader = angle_between(loads(p1), loads(p2), loads(p3))
                if sving_grader > 350 or sving_grader < 10:
                    veger_edges = veger_edges.drop(idx, axis=0)
                    
# %%
