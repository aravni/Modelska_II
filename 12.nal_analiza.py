# coding: utf-8
import numpy as np 
import matplotlib.pyplot as plt

Re,N_kv,iter,t,Fy_n,Fy_s,Fx_w,Fx_e,continuity,alpha_p,alpha_v = np.loadtxt('./results_SIMPLE.txt', skiprows=2, usecols=(0,2,3,4,5,6,7,8,9,10,11),unpack=True)

################# OK -done
# for i in range(len(iter)):
#     if Re[i] == 1:
#         plt.plot(iter[i],t[i],'C0o',label='1')
#     elif Re[i] == 10:
#         plt.plot(iter[i],t[i],'C1*',label='10')
#     elif Re[i] == 50:
#         plt.plot(iter[i],t[i],'C2>',label='50')
#     elif Re[i] == 100:
#         plt.plot(iter[i],t[i],'C3v',label='100')
#     elif Re[i] == 200:
#         plt.plot(iter[i],t[i],'C4s',label='200')
# legend = plt.legend(handles=[plt.Line2D([], [], color='C0', marker='o', label='1'),
#                              plt.Line2D([], [], color='C1', marker='*', label='10'),
#                              plt.Line2D([], [], color='C2', marker='>', label='50'),
#                              plt.Line2D([], [], color='C3', marker='v', label='100'),
#                              plt.Line2D([], [], color='C4', marker='s', label='200')],
#                              title='$Re$')
# plt.gca().add_artist(legend)
# plt.grid()
# plt.xlabel('$N_{iter}$')
# plt.ylabel('$t_{sim}$ [s]')
# plt.show()

################# OK - done
# for i in range(len(iter)):
#     if N_kv[i] == 2500:
#         if Re[i] == 1:
#             plt.plot(iter[i],continuity[i],'C0o',label='1')
#         elif Re[i] == 10:
#             plt.plot(iter[i],continuity[i],'C1*',label='10')
#         # elif Re[i] == 50:
#         #     plt.plot(iter[i],continuity[i],'C2>',label='50')
#         elif Re[i] == 100:
#             plt.plot(iter[i],continuity[i],'C3v',label='100')
#         elif Re[i] == 200:
#             plt.plot(iter[i],continuity[i],'C4s',label='200')
# legend = plt.legend(handles=[plt.Line2D([], [], color='C0', marker='o', label='1'),
#                              plt.Line2D([], [], color='C1', marker='*', label='10'),
#                              plt.Line2D([], [], color='C3', marker='v', label='100'),
#                              plt.Line2D([], [], color='C4', marker='s', label='200')],
#                              title='$Re$',loc='best')
# plt.gca().add_artist(legend)
# plt.grid()
# plt.yscale('log')
# plt.xlabel('$N_{iter}$')
# plt.ylabel('mass imbalance residual')
# plt.title('$N_{KV}=2500$')
# plt.show()

################# OK - done
# for i in range(len(iter)):
#     if Re[i] == 1:
#         plt.plot(N_kv[i],continuity[i],'C0o')
# #     elif Re[i] == 10:
# #         plt.plot(N_kv[i],continuity[i],'C1*')
# #     elif Re[i] == 50:
# #         plt.plot(N_kv[i],continuity[i],'C2>')
# #     elif Re[i] == 100:
# #         plt.plot(N_kv[i],continuity[i],'C3v')
# #     elif Re[i] == 200:
# #         plt.plot(N_kv[i],continuity[i],'C4s')
# # legend = plt.legend(handles=[plt.Line2D([], [], color='C0', marker='o', label='1'),
# #                              plt.Line2D([], [], color='C1', marker='*', label='10'),
# #                              plt.Line2D([], [], color='C2', marker='>', label='50'),
# #                              plt.Line2D([], [], color='C3', marker='v', label='100'),
# #                              plt.Line2D([], [], color='C4', marker='s', label='200')],
# #                              title='$Re$')
# # plt.gca().add_artist(legend)
# plt.grid()
# plt.yscale('log')
# plt.xlabel('$N_{KV}$')
# plt.ylabel('mass imbalance residual')
# plt.title('$Re=1$')
# plt.show()

################# OK - done
# for i in range(len(iter)):
#     if Re[i]>1 and iter[i]>600:
#         plt.plot(Re[i],Fy_n[i],'C0o')
#         plt.plot(Re[i],Fy_s[i],'C1o')
#         plt.plot(Re[i],abs(Fx_w[i]),'C2o')
#         plt.plot(Re[i],Fx_e[i],'C3o')
#     # elif iter[i] == 10:
#     #     plt.plot(Re[i],Fy_n[i],'C1*')
#     # elif iter[i] == 50:
#     #     plt.plot(Re[i],Fy_n[i],'C2>')
#     # elif iter[i] == 100:
#     #     plt.plot(Re[i],Fy_n[i],'C3v')
#     # elif iter[i] == 200:
#     #     plt.plot(Re[i],Fy_n[i],'C4s')
# legend = plt.legend(handles=[plt.Line2D([], [], color='C0', label='zgoraj'),
#                              plt.Line2D([], [], color='C1', label='spodaj'),
#                              plt.Line2D([], [], color='C2', label='levo'),
#                              plt.Line2D([], [], color='C3', label='desno')],
#                              title='rob kotanje')
# plt.gca().add_artist(legend)
# plt.grid()
# plt.yscale('log')
# plt.xlabel('$Re$')
# plt.ylabel('$F$')
# plt.show()

################# OK - done
# for i in range(len(iter)):
#     if Re[i] == 1 and N_kv[i] < 500:
#         plt.plot(iter[i],Fy_n[i],'C0o')
#     # elif Re[i] == 10:
#     #     plt.plot(iter[i],Fy_n[i],'C1*')
#     # elif Re[i] == 50:
#     #     plt.plot(iter[i],Fy_n[i],'C2>')
#     elif Re[i] == 100:
#         plt.plot(iter[i],Fy_n[i],'C3v')
#     # elif Re[i] == 200:
#     #     plt.plot(iter[i],Fy_n[i],'C4s')
# legend = plt.legend(handles=[plt.Line2D([], [], color='C0', marker='o', label='1'),
#                              plt.Line2D([], [], color='C3', marker='v', label='100')],
#                              title='$Re$')
# plt.gca().add_artist(legend)
# plt.grid()
# plt.xlabel('$N_{iter}$')
# plt.ylabel('sila na pokrov')
# plt.show()

################# OK - done
# for i in range(len(iter)):
#     if Re[i] == 1:
#         plt.plot(N_kv[i],Fy_n[i],'C0o')
#         plt.plot(N_kv[i],Fy_s[i],'C1*')
#         plt.plot(N_kv[i],abs(Fx_w[i]),'C2s')
#         plt.plot(N_kv[i],Fx_e[i],'C3v')
#     # elif Re[i] == 10:
#     #     plt.plot(N_kv[i],Fy_n[i],'C1*')
#     # elif Re[i] == 50:
#     #     plt.plot(N_kv[i],Fy_n[i],'C2>')
#     # elif Re[i] == 100:
#     #     plt.plot(N_kv[i],Fy_n[i],'C3v')
#     # elif Re[i] == 200:
#     #     plt.plot(N_kv[i],Fy_n[i],'C4s')
# legend = plt.legend(handles=[plt.Line2D([], [], color='C0', marker='o', label='zgoraj'),
#                              plt.Line2D([], [], color='C1', marker='*', label='spodaj'),
#                              plt.Line2D([], [], color='C2', marker='s', label='levo'),
#                              plt.Line2D([], [], color='C3', marker='v', label='desno')],
#                              title='rob kotanje')
# plt.gca().add_artist(legend)
# plt.grid()
# plt.yscale('log')
# plt.title('$Re=1$')
# plt.xlabel('$N_{kv}$')
# plt.ylabel('$F$')
# plt.show()