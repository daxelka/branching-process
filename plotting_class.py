import matplotlib.pyplot as plt


class Plotting:

    def blue_orange_green_set(self):
        return '#1f77b4', '#ff7f0e', '#2ca02c'

    def purple_red_teal_set(self):
        # return '#9467bd', '#8c564b', '#17becf'
        return '#9467bd', '#af3636', '#17becf', '#FFDB58', '#A9A9A9'

    def CB91_color_set(self):
        return '#2CBDFE', '#47DBCD', '#F3A0F2', '#9D2EC5', '#661D98', '#F5B14C'

    def plot(self, curves_list, title='my title', xlim=[0, 10], if_safe=False):
        my_colors = self.purple_red_teal_set()
        my_linewidth = 0.8
        my_fontsize = 14
        my_mark_size = 14

        for n, curves_dict in enumerate(curves_list):
            plt.plot(curves_dict['gens_analyt'], curves_dict['probs_analyt'],
                     c=my_colors[n], linewidth=my_linewidth, label=curves_dict['label_analyt'])
            plt.scatter(curves_dict['gens_num'], curves_dict['probs_num'],
                        c=my_colors[n], s=my_mark_size, label=curves_dict['label_num'])

        plt.xlabel('generation', fontsize=my_fontsize)
        plt.ylabel('probability', fontsize=my_fontsize)
        plt.title(title, fontsize=my_fontsize)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.legend(frameon=False)
        plt.xlim(xlim)
        if if_safe:
            plt.savefig('img/' + title + '.png', dpi=300)
        else:
            plt.show()
        print('plotted')