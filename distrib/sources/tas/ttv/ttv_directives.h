
void ttv_setdirectives (ttvfig_list *tvf, inffig_list *ifl, HeapAlloc *ha);
ttv_directive *ttv_get_directive(ttvsig_list *tvs);
chain_list *ttv_get_signal_with_directives(ttvfig_list *tvf, int clocks, int data);
