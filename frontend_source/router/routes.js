
const routes = [
  {
    path: '/',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', component: () => import('pages/MoleculePage.vue') },
      { path: '/reaction', component: () => import('pages/ReactionPage.vue') },
      { path: '/settings', component: () => import('pages/SettingsPage.vue') },
      { path: '/tutorial', component: () => import('pages/TutorialPage.vue') },
      { path: '/quintessen-demo', component: () => import('pages/QuintessenLabsDemo.vue') }
    ]
  },

  // Always leave this as last one,
  // but you can also remove it
  {
    path: '/:catchAll(.*)*',
    component: () => import('pages/ErrorNotFound.vue')
  }
]

export default routes
