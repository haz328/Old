from ESP_Class import *


# parametric study for surging test
def surging_test():
    conn, c = connect_db('ESP.db')

    te2700_case = TwoPhaseCompare('TE2700', conn)
    df_50psi = te2700_case.surging_performance(2700, 0.2, 3500, 50, 100)
    df_100psi = te2700_case.surging_performance(2700, 0.2, 3500, 100, 100)
    df_150psi = te2700_case.surging_performance(2700, 0.2, 3500, 150, 100)
    df_200psi = te2700_case.surging_performance(2700, 0.2, 3500, 200, 100)

    dfp = pd.DataFrame({'50psi': df_50psi.zhu,
                        '100psi': df_100psi.zhu,
                        '150psi': df_150psi.zhu,
                        '200psi': df_200psi.zhu,
                        'gf_50psi': df_50psi.gf,
                        'gf_100psi': df_100psi.gf,
                        'gf_150psi': df_150psi.gf,
                        'gf_200psi': df_200psi.gf})

    df_3500rpm = te2700_case.surging_performance(2700, 0.2, 3500, 150, 100)
    df_3000rpm = te2700_case.surging_performance(2700 * 3/3.5, 0.2, 3000, 150, 100)
    df_2500rpm = te2700_case.surging_performance(2700 * 2.5/3.5, 0.2, 2500, 150, 100)
    df_2000rpm = te2700_case.surging_performance(2700 * 2/3.5, 0.2, 2000, 150, 100)

    dfn = pd.DataFrame({'3500rpm': df_3500rpm.zhu,
                        '3000rpm': df_3000rpm.zhu,
                        '2500rpm': df_2500rpm.zhu,
                        '2000rpm': df_2000rpm.zhu,
                        'gf_3500rpm': df_3500rpm.gf,
                        'gf_3000rpm': df_3000rpm.gf,
                        'gf_2500rpm': df_2500rpm.gf,
                        'gf_2000rpm': df_2000rpm.gf})

    fig1 = plt.figure(dpi=300)
    ax1 = fig1.add_subplot(111)
    ax1.plot(dfp['gf_50psi'] * 100, dfp['50psi'], 'mo', label='50 psig')
    ax1.plot(dfp['gf_100psi'] * 100, dfp['100psi'], 'ks', label='100 psig')
    ax1.plot(dfp['gf_150psi'] * 100, dfp['150psi'], 'r^', label='150 psig')
    ax1.plot(dfp['gf_200psi'] * 100, dfp['200psi'], 'g*', label='200 psig')

    ax1.set_xlabel(r'$\lambda_{G}$ (%)')
    ax1.set_ylabel(r'$N_{P}$')
    ax1.legend(frameon=False)

    fig2 = plt.figure(dpi=300)
    ax2 = fig2.add_subplot(111)
    ax2.plot(dfn['gf_3500rpm'] * 100, dfn['3500rpm'], 'mo', label='3500 rpm')
    ax2.plot(dfn['gf_3000rpm'] * 100, dfn['3000rpm'], 'ks', label='3000 rpm')
    ax2.plot(dfn['gf_2500rpm'] * 100, dfn['2500rpm'], 'r^', label='2500 rpm')
    ax2.plot(dfn['gf_2000rpm'] * 100, dfn['2000rpm'], 'g*', label='2000 rpm')

    ax2.set_xlabel(r'$\lambda_{G}$ (%)')
    ax2.set_ylabel(r'$N_{P}$')
    ax2.legend(frameon=False)

    fig1.show()
    fig2.show()
    disconnect_db(conn)


# parametric study for mapping test
def mapping_test():
    conn, c = connect_db('ESP.db')

    te2700_case = TwoPhaseCompare('TE2700', conn)
    df30psi = te2700_case.mapping_performance(100, 4600, 3500, 30, 100)
    df100psi = te2700_case.mapping_performance(100, 4600, 3500, 100, 100)
    df500psi = te2700_case.mapping_performance(100, 4600, 3500, 500, 100)
    df1000psi = te2700_case.mapping_performance(100, 4600, 3500, 1000, 100)

    dfp = pd.DataFrame({'30psi': df30psi.zhu,
                        '100psi': df100psi.zhu,
                        '500psi': df500psi.zhu,
                        '1000psi': df1000psi.zhu,
                        'ql_30psi': df30psi.ql,
                        'ql_100psi': df100psi.ql,
                        'ql_500psi': df500psi.ql,
                        'ql_1000psi': df1000psi.ql})

    df3500rpm = te2700_case.mapping_performance(100, 4600, 3500, 150, 100)
    df3000rpm = te2700_case.mapping_performance(100, 4600 * 3/3.5, 3000, 150, 100)
    df2500rpm = te2700_case.mapping_performance(100, 4600 * 2.5/3.5, 2500, 150, 100)
    df2000rpm = te2700_case.mapping_performance(100, 4600 * 2/3.5, 2000, 150, 100)

    dfn = pd.DataFrame({'3500rpm': df3500rpm.zhu,
                        '3000rpm': df3000rpm.zhu,
                        '2500rpm': df2500rpm.zhu,
                        '2000rpm': df2000rpm.zhu,
                        'ql_3500rpm': df3500rpm.ql,
                        'ql_3000rpm': df3000rpm.ql,
                        'ql_2500rpm': df2500rpm.ql,
                        'ql_2000rpm': df2000rpm.ql})

    fig1 = plt.figure(dpi=300)
    ax1 = fig1.add_subplot(111)
    ax1.plot(dfp['ql_30psi'], dfp['30psi'], 'g*', label='30 psig')
    ax1.plot(dfp['ql_100psi'], dfp['100psi'], 'r^', label='100 psig')
    ax1.plot(dfp['ql_500psi'], dfp['500psi'], 'ks', label='500 psig')
    ax1.plot(dfp['ql_1000psi'], dfp['1000psi'], 'mo', label='1000 psig')

    ax1.set_xlabel(r'$Q_{L}$ (bpd)')
    ax1.set_ylabel(r'$P$ (psi)')
    ax1.set_xlim(0, 5000)
    ax1.legend(frameon=False)

    fig2 = plt.figure(dpi=300)
    ax2 = fig2.add_subplot(111)
    ax2.plot(dfn['ql_3500rpm'], dfn['3500rpm'], 'g*', label='3500rpm')
    ax2.plot(dfn['ql_3000rpm'], dfn['3000rpm'], 'r^', label='3000rpm')
    ax2.plot(dfn['ql_2500rpm'], dfn['2500rpm'], 'ks', label='2500rpm')
    ax2.plot(dfn['ql_2000rpm'], dfn['2000rpm'], 'mo', label='2000rpm')

    ax2.set_xlabel(r'$Q_{L}$ (bpd)')
    ax2.set_ylabel(r'$P$ (psi)')
    ax2.set_xlim(0, 5000)
    ax2.legend(frameon=False)

    fig1.show()
    fig2.show()
    disconnect_db(conn)


if __name__ == "__main__":
    # surging_test()
    mapping_test()

