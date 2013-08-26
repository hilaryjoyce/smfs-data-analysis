# Visualizing clusters code

def cluster_tf_list(curve_list, coa):
    from essentials import load2Col
    tf_files = [coa.tss_force_files[x] for x in curve_list]
    tf_list = []
    for file in tf_files:
        tf_list.append(load2Col(file))
    return tf_list

def cluster_lcd_list(curve_list, coa):
    from essentials import load2Col
    coa_files = coa.list_density_files()
    lcd_files = [coa_files[x] for x in curve_list]
    lcd_list = []
    for file in lcd_files:
        lcd_list.append(load2Col(file))
    return lcd_list

def plot_subplots(data_list, id_list, max_x):
    from numpy import ceil
    import matplotlib.pyplot as plt
    cols = 4
    rows = ceil(len(data_list)/4.0)
    height = rows * 2.2
    width = cols * 3
    i = 1
    plt.figsize(width, height)
    while i <= len(data_list):
        data = data_list[i-1];
        plt.subplot(rows, cols, i);
        plt.plot(data[0], data[1]);
        plt.xlim(0, max_x);
        i = i+1;
    plt.tight_layout(pad=2);
    plt.subplots_adjust(top=0.96);
    return plt.gcf();

def plot_tf_subplots(sub_list, coa, max_x=100, max_y=300):
    from numpy import ceil
    import matplotlib.pyplot as plt
    data_list = cluster_tf_list(sub_list, coa)
    cols = 4
    rows = ceil(len(data_list)/4.0)
    height = rows * 2.2
    width = cols * 3
    i = 1
    plt.figsize(width, height)
    while i <= len(data_list):
        data = data_list[i-1];
        plt.subplot(rows, cols, i);
        plt.plot(data[0], data[1]);
        plt.xlim(0, max_x);
        plt.ylim(-200,max_y)
        plt.text(max_x/1.2, max_y/1.2, '%d' % i)
        i = i+1;
    plt.tight_layout(pad=2);
    plt.subplots_adjust(top=0.96);
    return plt.gcf();

def plot_lcd_subplots(data_list, coa, max_x=100, max_y=0.2):
    from numpy import ceil
    import matplotlib.pyplot as plt
    data_list = cluster_lcd_list(data_list, coa)
    cols = 4
    rows = ceil(len(data_list)/4.0)
    height = rows * 2.2
    width = cols * 3
    i = 1
    plt.figsize(width, height)
    while i <= len(data_list):
        data = data_list[i-1];
        plt.subplot(rows, cols, i);
        plt.plot(data[0], data[1]);
        plt.xlim(0, max_x);
        plt.ylim(0, max_y);
        plt.text(max_x/1.2, max_y/1.2, '%d' % i)
        i = i+1;
    plt.tight_layout(pad=2);
    plt.subplots_adjust(top=0.96);
    return plt.gcf();

def plot_lcd_tf_subplots(sub_list, coa, max_x, lcd_y=0.15, tf_y = 300, min_lc= 0):
    from numpy import ceil
    import matplotlib.pyplot as plt
    cpool = ["#1F78B4", "#E31A1C"]
    lcd_list = cluster_lcd_list(sub_list, coa)
    tf_list = cluster_tf_list(sub_list, coa)
    cols = 3
    rows = ceil(len(lcd_list)/3.0)
    height = rows * 2.2
    width = cols * 4
    i = 1
    plt.figsize(width, height)
    while i <= len(lcd_list):
        lcd = lcd_list[i-1];
        tf = tf_list[i-1];
        lc_ax = plt.subplot(rows, cols, i);
        lc_ax.plot(lcd[0], lcd[1], color=cpool[1], lw=1);
        lc_ax.set_ylim(min_lc, lcd_y);
        lc_ax.text(max_x*0.9, lcd_y*0.85, '%d' % i, va='top', ha='right')
        tf_ax = lc_ax.twinx()
        tf_ax.plot(tf[0], tf[1], color=cpool[0], alpha=0.6)
        tf_ax.set_ylim(-100, tf_y)
        lc_ax.set_xlim(0, max_x);
        tf_ax.set_yticks([])
        i = i+1;
        
    plt.tight_layout(pad=2);
    plt.subplots_adjust(top=0.96);
    return plt.gcf();
